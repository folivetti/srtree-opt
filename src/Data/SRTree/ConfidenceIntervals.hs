{-# language ViewPatterns #-}
module Data.SRTree.ConfidenceIntervals where

import qualified Numeric.LinearAlgebra as LA
import Statistics.Distribution hiding (Distribution)
import Statistics.Distribution.StudentT
import qualified Data.Vector.Storable as VS
import Data.SRTree
import Data.SRTree.Recursion ( cata )
import Data.SRTree.AD
import Data.SRTree.Likelihoods

data CIType = Laplace | Profile deriving (Eq, Show, Read)

data BasicStats = MkStats { _cov    :: LA.Matrix Double
                          , _corr   :: LA.Matrix Double
                          , _stdErr :: LA.Vector Double
                          } deriving (Eq, Show, Read)

data CI = CI { est_   :: Double
             , lower_ :: Double
             , upper_ :: Double
             } deriving (Eq, Show, Read)

data ProfileT = ProfileT { _taus   :: LA.Vector Double
                         , _thetas :: LA.Vector Double
                         , _deltas :: LA.Vector Double
                         } deriving (Eq, Show, Read)

getProfile dist mSErr xss ys tree theta ix stdErr = ProfileT taus thetas deltas
  where
    delta                      = stdErr_i / 8
    stdErr_i                   = stdErr LA.! ix
    theta_i                    = theta VS.! ix
    (taus, thetas, deltas)     = extractProfile $ go mempty 30 (-delta) 0 1 <> go mempty 30 delta
    nll_opt                    = nll dist mSErr xss ys tree theta
    tau_max = 10

    go 0 _     _ _         acc = acc
    go k delta t inv_slope acc
      | nll_cond < nll_opt = error "found a better theta"
      | abs tau > tau_max  = acc
      | otherwise          =  go ((tau, theta_t, delta_t) : acc) (k-1) stp t' inv_slope'
      where
        t'         = t + inv_slope
        inv_slope' = adjust $ abs (tau / (stdErr_i * zv))
        adjust x   = min 4.0 (max 0.0625 x)
        optTree    = replaceParamWithConst (theta_i + delta * t) ix tree
        theta'     = deleteAt ix theta
        theta_t    = insertAt ix (theta_i + delta * t) 
                   $ minimizeNLL dist mSErr 10 xss ys optTree theta'
        nll_cond   = nll dist mSErr xss ys tree theta_t
        delta_t    = (nll_opt - nll_cond) / stdErr_i

getStatsFromModel :: Distribution -> Maybe Double -> Columns -> Column -> Fix SRTree -> VS.Vector Double -> BasicStats
getStatsFromModel dist mSErr xss ys tree theta = MkStats cov corr stdErr
  where
    cov    = LA.inv $ LA.fromLists $ hessianNLL dist mSErr xss ys tree theta
    stdErr = sqrt $ LA.takeDiag cov
    corr   = cov / LA.outer stdErr stdErr

paramCI :: CIType -> Int -> LA.Vector Double -> VS.Vector Double -> Double -> [CI]
paramCI Laplace nSamples stdErr theta alpha = zipWith3 CI (VS.toList theta) lows highs
  where
    k     = VS.length theta
    t     = quantile (studentT . fromIntegral $ nSamples - k) (alpha / 2.0)
    lows  = VS.toList $ VS.zipWith (-) theta $ VS.map (*t) stdErr
    highs = VS.toList $ VS.zipWith (+) theta $ VS.map (*t) stdErr

predictionCI :: CIType -> LA.Matrix Double -> Columns ->  Fix SRTree -> VS.Vector Double -> Double -> [CI]
predictionCI Laplace cov xss tree theta alpha = zipWith3 CI (VS.toList yhat) lows highs
  where
    (yhat, LA.toRows . LA.fromColumns -> grad) = reverseModeUnique xss theta VS.singleton tree
    t = quantile (studentT . fromIntegral $ LA.size yhat - k) (alpha / 2.0)
    k = VS.length theta
    lows = VS.toList $ VS.zipWith (-) yhat $ VS.map (*t) resStdErr
    highs = VS.toList $ VS.zipWith (+) yhat $ VS.map (*t) resStdErr
    getResStdError row = sqrt $ LA.dot row $ LA.fromList $ map (LA.dot row . (cov LA.!)) [0 .. k-1]
    resStdErr          = VS.fromList $ map getResStdError grad

-- * Utilities (to be moved to SRTree)
replaceParamWithConst :: Int -> Double -> Fix SRTree -> Fix SRTree
replaceParamWithConst ix val = cata alg
  where
    alg (Param iy)   = if ix == iy then Const val else Param iy
    alg (Const c)    = Const c
    alg (Var iy)     = Var iy
    alg (Uni f t)    = Uni f t
    alg (Bin op l r) = Bin op l r
{-# inline replaceParamWithConst #-}

insertAt :: Int -> a -> VS.Vector a -> VS.Vector a
insertAt ix val = VS.fromList . go 0 . VS.toList
  where
      go iy [] = [val]
      go iy (y:ys)
        | iy == ix  = val : y : ys
        | otherwise = y : go (iy+1) ys
{-# inline insertAt #-}

deleteAt :: Int -> VS.Vector a -> VS.Vector a
deleteAt ix = VS.fromList . go 0 . VS.toList
  where
      go iy [] = []
      go iy (y:ys)
        | iy == ix  = y : ys
        | otherwise = y : go (iy+1) ys
{-# inline deleteAt #-}
