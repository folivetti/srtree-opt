{-# language ViewPatterns, ScopedTypeVariables #-}
module Data.SRTree.ConfidenceIntervals where

import qualified Numeric.LinearAlgebra as LA
import Statistics.Distribution hiding (Distribution)
import Statistics.Distribution.StudentT
import Statistics.Distribution.FDistribution
import qualified Data.Vector.Storable as VS
import qualified Data.Vector as V
import Data.SRTree
import Data.SRTree.Eval
import Data.SRTree.Recursion ( cata )
import Data.SRTree.AD
import Data.SRTree.Likelihoods
import Data.SRTree.Opt
import Numeric.GSL.Interpolation
import Data.List ( sortOn )
import qualified Numeric.LinearAlgebra.Static as LAS
import GHC.TypeLits
import Data.Proxy

import Debug.Trace ( trace )

data CIType = Laplace BasicStats | Profile BasicStats [ProfileT]

data BasicStats = MkStats { _cov    :: LA.Matrix Double
                          , _corr   :: LA.Matrix Double
                          , _stdErr :: LA.Vector Double
                          } deriving (Eq, Show, Read)

data CI = CI { est_   :: Double
             , lower_ :: Double
             , upper_ :: Double
             } deriving (Eq, Show, Read)

data ProfileT = ProfileT { _taus      :: LA.Vector Double
                         , _thetas    :: LA.Matrix Double
                         , _deltas    :: LA.Vector Double
                         , _tau2theta :: Double -> Double
                         , _theta2tau :: Double -> Double
                         }

showCI :: Int -> CI -> String
showCI n (CI x l h) = show (rnd l) <> " <= " <> show (rnd x) <> " <= " <> show (rnd h)
  where
      rnd = (/10^n) . (fromIntegral . round) . (*10^n)
printCI :: Int -> CI -> IO ()
printCI n = putStrLn . showCI n

paramCI :: CIType -> Int -> VS.Vector Double -> Double -> [CI]
paramCI (Laplace stats) nSamples theta alpha = zipWith3 CI (VS.toList theta) lows highs
  where
    k      = VS.length theta
    t      = quantile (studentT . fromIntegral $ nSamples - k) (1 - alpha / 2.0)
    stdErr = _stdErr stats
    lows   = VS.toList $ VS.zipWith (-) theta $ VS.map (*t) stdErr
    highs  = VS.toList $ VS.zipWith (+) theta $ VS.map (*t) stdErr

paramCI (Profile stats profiles) nSamples theta alpha = zipWith3 CI (VS.toList theta) lows highs
  where
    k        = VS.length theta
    t        = quantile (studentT . fromIntegral $ nSamples - k) (1 - alpha / 2.0)
    stdErr   = _stdErr stats
    lows     = map (\p -> _tau2theta p (-t)) profiles
    highs    = map (\p -> _tau2theta p t) profiles

-- predictionCI :: CIType -> Columns -> Column -> Columns -> Fix SRTree -> VS.Vector Double -> Double -> [CI]
predictionCI :: CIType -> (Columns -> Column) -> (Columns -> [VS.Vector Double]) -> (VS.Vector Double -> Fix SRTree -> Double -> Double) -> Columns -> Fix SRTree -> VS.Vector Double -> Double -> [CI]
predictionCI (Laplace stats) predFun jacFun _ xss tree theta alpha = zipWith3 CI (VS.toList yhat) lows highs
  where
    yhat  = predFun xss
    jac   = jacFun xss
    t     = quantile (studentT . fromIntegral $ LA.size yhat - k) (1 - alpha / 2.0)
    cov   = _cov stats
    k     = VS.length theta
    lows  = VS.toList $ VS.zipWith (-) yhat $ VS.map (*t) resStdErr
    highs = VS.toList $ VS.zipWith (+) yhat $ VS.map (*t) resStdErr

    getResStdError row = sqrt $ LA.dot row $ LA.fromList $ map (LA.dot row . (cov LA.!)) [0 .. k-1]
    resStdErr          = VS.fromList $ map getResStdError jac

predictionCI (Profile _ _) predFun _ profFun xss tree theta alpha = zipWith f (VS.toList yhat) xss'
  where
    yhat = predFun xss

    t       = quantile (studentT . fromIntegral $ LA.size yhat - VS.length theta) (alpha / 2.0)
    theta0  = calcTheta0 tree
    xss'    = LA.toRows $ LA.fromColumns $ V.toList xss
    
    f yh xs = let t' = replaceParam0 tree $ evalVar xs theta0
                  spline = profFun theta t'
               in CI yh (spline t) (spline (-t))

replaceParam0 :: Fix SRTree -> Fix SRTree -> Fix SRTree
replaceParam0 tree t0 = cata alg tree
  where
    alg (Var ix)     = Fix $ Var ix
    alg (Param 0)    = t0
    alg (Param ix)   = Fix $ Param ix
    alg (Const c)    = Fix $ Const c
    alg (Uni g t)    = Fix $ Uni g t
    alg (Bin op l r) = Fix $ Bin op l r

evalVar :: VS.Vector Double -> Fix SRTree -> Fix SRTree
evalVar xs = cata alg 
  where
    alg (Var ix)     = Fix $ Const (xs VS.! ix)
    alg (Param ix)   = Fix $ Param ix
    alg (Const c)    = Fix $ Const c
    alg (Uni g t)    = Fix $ Uni g t
    alg (Bin op l r) = Fix $ Bin op l r

calcTheta0 :: Fix SRTree -> Fix SRTree
calcTheta0 tree = case cata alg tree of
                         Left g -> g (Fix $ Param 0)
                         Right v -> error "No theta0?"
  where
    alg (Var ix)     = Right $ Fix $ Var ix
    alg (Param 0)    = Left id
    alg (Param ix)   = Right $ Fix $ Param ix
    alg (Const c)    = Right $ Fix $ Const c
    alg (Uni g t)    = case t of
                         Left f  -> Left $ f . evalInverse g
                         Right v -> Right $ evalFun g v
    alg (Bin op l r) = case l of
                         Left f   -> case r of
                                       Left  g -> error "This shouldn't happen!"
                                       Right v -> Left $ f . invright op v
                         Right vl -> case r of
                                       Left  g -> Left $ g . invleft op vl
                                       Right vr -> Right $ evalOp op vl vr

getAllProfiles dist mSErr xss ys tree theta stdErr = reverse (getAll 0 [])
  where 
    nParams    = VS.length theta
    tau_max    = sqrt $ quantile (fDistribution nParams (LA.size ys - nParams)) (1 - 0.01)
    profFun ix = getProfile dist mSErr xss ys tree theta (stdErr LA.! ix) tau_max ix

    getAll ix acc | ix == nParams = acc
                  | otherwise     = case profFun ix of
                                      Left t  -> getAllProfiles dist mSErr xss ys tree t stdErr
                                      Right p -> getAll (ix + 1) (p : acc)

getProfile dist mSErr xss ys tree theta stdErr_i tau_max ix = 
  do negDelta <- go 30 (-stdErr_i / 8) 0 1 mempty
     posDelta <- go 30  (stdErr_i / 8) 0 1 ([0], [theta], [0])
     let (LA.fromList -> taus, LA.fromRows -> thetas, LA.fromList -> deltas) = negDelta <> posDelta
         (tau2theta, theta2tau) = createSplines taus thetas stdErr_i ix
     pure $ ProfileT taus thetas deltas tau2theta theta2tau
  where
    nll_opt = nll dist mSErr xss ys tree theta
    go 0 _     _ _         acc = Right acc
    go k delta t inv_slope acc
      | nll_cond < nll_opt = Left theta_t
      | abs tau > tau_max  = Right (([tau], [theta_t], [delta_t]) <> acc) 
      | otherwise          = go (k-1) delta (t + inv_slope) inv_slope' (([tau], [theta_t], [delta_t]) <> acc) 
      where
        zv         = gradNLL dist mSErr xss ys tree theta_t VS.! ix
        inv_slope' = min 4.0 . max 0.0625 . abs $ (tau / (stdErr_i * zv))
        theta_t    = fst
                   $ minimizeNLLWithFixedParam dist mSErr 10 xss ys tree ix 
                   $ theta VS.// [(ix, (theta VS.! ix) + delta * (t + inv_slope))]
        nll_cond   = nll dist mSErr xss ys tree theta_t
        tau        = signum delta * sqrt (2*nll_cond - 2*nll_opt)
        delta_t    = (nll_cond - nll_opt) / stdErr_i

getStatsFromModel :: Distribution -> Maybe Double -> Columns -> Column -> Fix SRTree -> VS.Vector Double -> BasicStats
getStatsFromModel dist mSErr xss ys tree theta = MkStats cov corr stdErr
  where
    nParams = fromIntegral $ LA.size theta
    hess :: LA.Matrix Double
    hess    = LA.fromLists $ hessianNLL dist mSErr xss ys tree theta
    cov     = case LA.mbChol (LA.sym hess) of
                Nothing -> error "Hessian of the model is not positive definite"
                Just m  -> LA.cholSolve m (LA.ident nParams)
    
    stdErr = sqrt $ LA.takeDiag cov
    corr   = cov / LA.outer stdErr stdErr

-- Create splines for profile-t
createSplines :: VS.Vector Double -> LA.Matrix Double -> Double -> Int -> (Double -> Double, Double -> Double)
createSplines taus thetas se ix
  | VS.length taus < 2 = (evaluate CSpline [(0, -se), (1, se)], evaluate CSpline [(-se, 0), (se, 1)])
  | otherwise          = (tau2theta, theta2tau)
  where
    tau2theta = evaluate CSpline $ sortOnFirst taus (getCol ix thetas)
    theta2tau = evaluate CSpline $ sortOnFirst (getCol ix thetas) taus

getCol :: Int -> LA.Matrix Double -> LA.Vector Double
getCol ix mtx = LA.flatten $ mtx LA.¿ [ix]
{-# inline getCol #-}

sortOnFirst :: VS.Vector Double -> VS.Vector Double -> [(Double, Double)]
sortOnFirst xs ys = sortOn fst $ zip (VS.toList xs) (VS.toList ys)
{-# inline sortOnFirst #-}

splinesSketches :: Double -> VS.Vector Double -> VS.Vector Double -> (Double -> Double) -> (Double -> Double)
splinesSketches tauScale (VS.toList -> tau) (VS.toList -> theta) theta2tau
  | length tau < 2 = id
  | otherwise      = evaluate CSpline gpq
  where
    gpq = sortOn fst [(x, acos y') | (x, y) <- zip tau theta
                                   , let y' = theta2tau y / tauScale
                                   , abs y' < 1 ]

approximateContour :: Int -> Int -> [ProfileT] -> Int -> Int -> Double -> [(Double, Double)]
approximateContour nParams nPoints profs ix1 ix2 alpha = go 0
  where
    -- get the info for ix1 and ix2
    (prof1, prof2)           = (profs !! ix1, profs !! ix2)
    (tau2theta1, theta2tau1) = (_tau2theta prof1, _theta2tau prof1)
    (tau2theta2, theta2tau2) = (_tau2theta prof2, _theta2tau prof2)

    -- calculate the spline for A-D
    tauScale = sqrt (fromIntegral nParams * quantile (fDistribution nParams (nPoints - nParams)) (1 - alpha))
    splineG1 = splinesSketches tauScale (_taus prof1) (getCol ix2 (_thetas prof1)) theta2tau2
    splineG2 = splinesSketches tauScale (_taus prof2) (getCol ix1 (_thetas prof2)) theta2tau1
    angles   = [ (0, splineG1 1), (splineG2 1, 0), (pi, splineG1 (-1)), (splineG2 (-1), pi) ]
    splineAD = evaluate CSpline points

    applyIfNeg (x, y) = if y < 0 then (-x, -y) else (x ,y)
    points   = sortOn fst
             $ [applyIfNeg ((x+y)/2, x - y) | (x, y) <- angles] 
            <> (\(x,y) -> [(x + 2*pi, y)]) (head points)

    -- generate the points of the curve
    go 100 = []
    go ix  = (p, q) : go (ix+1)
      where
        ai = ix * 2 * pi / 99 - pi
        di = splineAD ai
        taup = cos (ai + di / 2) * tauScale
        tauq = cos (ai - di / 2) * tauScale
        p = tau2theta1 taup
        q = tau2theta2 tauq
