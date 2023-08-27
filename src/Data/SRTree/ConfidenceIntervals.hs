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

paramCI :: CIType -> Int -> VS.Vector Double -> Double -> [CI]
paramCI (Laplace stats) nSamples theta alpha = zipWith3 CI (VS.toList theta) lows highs
  where
    k      = VS.length theta
    t      = quantile (studentT . fromIntegral $ nSamples - k) (alpha / 2.0)
    stdErr = _stdErr stats
    lows   = VS.toList $ VS.zipWith (-) theta $ VS.map (*t) stdErr
    highs  = VS.toList $ VS.zipWith (+) theta $ VS.map (*t) stdErr

paramCI (Profile stats profiles) nSamples theta alpha = zipWith3 CI (VS.toList theta) lows highs
  where
    k        = VS.length theta
    t        = quantile (studentT . fromIntegral $ nSamples - k) (alpha / 2.0)
    stdErr   = _stdErr stats
    lows     = map (\p -> _tau2theta p (-t)) profiles
    highs    = map (\p -> _tau2theta p t) profiles

-- predictionCI :: CIType -> Columns -> Column -> Columns -> Fix SRTree -> VS.Vector Double -> Double -> [CI]
predictionCI (Laplace stats) predFun jacFun _ xss tree theta alpha = zipWith3 CI (VS.toList yhat) lows highs
  where
    yhat  = predFun xss
    jac   = jacFun xss
    t     = quantile (studentT . fromIntegral $ LA.size yhat - k) (alpha / 2.0)
    cov   = _cov stats
    k     = VS.length theta
    lows  = VS.toList $ VS.zipWith (-) yhat $ VS.map (*t) resStdErr
    highs = VS.toList $ VS.zipWith (+) yhat $ VS.map (*t) resStdErr

    getResStdError row = sqrt $ LA.dot row $ LA.fromList $ map (LA.dot row . (cov LA.!)) [0 .. k-1]
    resStdErr          = VS.fromList $ map getResStdError jac

{-
predFun xss = predict tree theta xss
jacFun  xss = LA.toRows . LA.fromColumns $ reverseModeUnique xss theta VS.singleton tree
profFun theta = let stdErr_0 = (VS.! 0) $  _stdErr $ getStatsFromModel dist mSErr xss_tr ys_tr tree theta
                 in _tau2theta $ getProfile dist mSErr xss_tr ys_tr tree theta stdErr_0 tau_max 0
-}
predictionCI (Profile _ _) predFun _ profFun xss tree theta alpha = zipWith f (VS.toList yhat) (VS.toList thetas0)
  where
    yhat = predFun xss

    t       = quantile (studentT . fromIntegral $ LA.size yhat - VS.length theta) (alpha / 2.0)
    thetas0 = calcTheta0 xss yhat theta VS.singleton tree
    tau_max = sqrt $ quantile (fDistribution (VS.length theta) (LA.size yhat - VS.length theta)) (1 - 0.01)
    
    f yh t0 = let spline = profFun (theta VS.// [(0, t0)])
               in CI yh (spline (-t)) (spline t)

calcTheta0 xss ys theta f tree = case cata alg tree of
                                              Left g -> g ys
                                              Right v -> error "No theta0?"
  where
    alg (Var ix)     = Right $ xss V.! ix
    alg (Param 0)    = Left id
    alg (Param ix)   = Right $ f $ theta VS.! ix
    alg (Const c)    = Right $ f c
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
      | abs tau > tau_max  = Right acc
      | otherwise          = go (k-1) delta (t + inv_slope) inv_slope' (([tau], [theta_t], [delta_t]) <> acc) 
      where
        zv         = gradNLL dist mSErr xss ys tree theta_t VS.! ix
        inv_slope' = min 4.0 . max 0.0625 . abs $ (tau / (stdErr_i * zv))
        theta_t    = fst
                   $ minimizeNLLWithFixedParam dist mSErr 10 xss ys tree ix 
                   $ theta VS.// [(ix, (theta VS.! ix) + delta * t)]
        nll_cond   = nll dist mSErr xss ys tree theta_t
        tau        = signum delta_t * sqrt (2*nll_cond - 2*nll_opt)
        delta_t    = (nll_opt - nll_cond) / stdErr_i

getStatsFromModel :: Distribution -> Maybe Double -> Columns -> Column -> Fix SRTree -> VS.Vector Double -> BasicStats
getStatsFromModel dist mSErr xss ys tree theta = MkStats cov corr stdErr
  where
    nParams = fromIntegral $ LA.size theta
    hess    = LA.fromLists $ hessianNLL dist mSErr xss ys tree theta
    {--
    hs      = case someNatVal nParams of
               Just (SomeNat (_ :: Proxy n)) -> case (LAS.create hess) :: Maybe (LAS.L n n) of
                 Nothing -> error "incorrect dimensions"
                 Just m  -> m :: LAS.L n n
               Nothing -> error "wat?"
    --}
    cov     = LA.inv hess -- hs LAS.<> (LAS.tr hs) -- LA.cholSolve (hs LA.<> (LA.tr hs)) (LA.ident nParams) -- LA.inv $ LA.fromLists $ hessianNLL dist mSErr xss ys tree theta
    stdErr = sqrt $ LA.takeDiag cov
    corr   = cov / LA.outer stdErr stdErr

-- Create splines for profile-t
createSplines :: VS.Vector Double -> LA.Matrix Double -> Double -> Int -> (Double -> Double, Double -> Double)
createSplines taus thetas se ix
  | VS.length taus < 2 = (evaluate CSpline [(0, -se), (1, se)], evaluate CSpline [(-se, 0), (se, 1)])
  | otherwise          = (tau2theta, theta2tau)
  where
    tau2theta = evaluate CSpline $ sortOnFirst taus (getRow ix thetas)
    theta2tau = evaluate CSpline $ sortOnFirst (getRow ix thetas) taus

getRow :: Int -> LA.Matrix Double -> LA.Vector Double
getRow ix mtx = LA.flatten $ mtx LA.? [ix]
{-# inline getRow #-}

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
    splineG1 = splinesSketches tauScale (_taus prof1) (getRow ix2 (_thetas prof1)) theta2tau2
    splineG2 = splinesSketches tauScale (_taus prof2) (getRow ix1 (_thetas prof2)) theta2tau1
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