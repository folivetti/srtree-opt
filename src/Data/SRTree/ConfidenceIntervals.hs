{-# language ViewPatterns #-}
module Data.SRTree.ConfidenceIntervals where

import qualified Numeric.LinearAlgebra as LA
import Statistics.Distribution hiding (Distribution)
import Statistics.Distribution.StudentT
import Statistics.Distribution.FDistribution
import qualified Data.Vector.Storable as VS
import Data.SRTree
import Data.SRTree.Recursion ( cata )
import Data.SRTree.AD
import Data.SRTree.Likelihoods
import Data.SRTree.Opt
import Numeric.GSL.Interpolation
import Data.List ( sortOn )

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
                         , _thetas :: LA.Matrix Double
                         , _deltas :: LA.Vector Double
                         } deriving (Eq, Show, Read)

getAllProfiles dist mSErr xss ys tree theta stdErr = go 0
  where 
    nParams = VS.length theta
    tau_max = sqrt $ quantile (fDistribution nParams (LA.size ys - nParams)) (1 - 0.01)
    go ix | ix == nParams = []
          | otherwise     = case getProfile dist mSErr xss ys tree theta stdErr tau_max ix of
                              Left t  -> getAllProfiles dist mSErr xss ys tree t stdErr
                              Right p -> p : go (ix + 1)

getProfile dist mSErr xss ys tree theta stdErr tau_max ix = 
  case go 30 (-delta) 0 1 mempty of
    Left t -> Left t
    Right negDelta -> case go 30 (-delta) 0 1 mempty of
                        Left t -> Left t
                        Right posDelta -> let (taus, thetas, deltas) = negDelta <> ([0], [theta], [0]) <> posDelta
                                           in Right $ ProfileT (LA.fromList taus) (LA.fromRows thetas) (LA.fromList deltas)
  where
    delta                      = stdErr_i / 8
    stdErr_i                   = stdErr LA.! ix

    theta_i                    = theta VS.! ix
    thetaCond                  = VS.take ix theta <> VS.tail (VS.drop ix theta)

    nll_opt                    = nll dist mSErr xss ys tree theta


    go :: Int -> Double -> Double -> Double -> ([Double], [VS.Vector Double], [Double]) -> Either (VS.Vector Double) ([Double], [VS.Vector Double], [Double])
    go 0 _     _ _         acc = Right acc
    go k delta t inv_slope acc
      | nll_cond < nll_opt = Left theta_t
      | abs tau > tau_max  = Right acc
      | otherwise          = go (k-1) delta t' inv_slope' (([tau], [theta_t], [delta_t]) <> acc) 
      where
        t'         = t + inv_slope
        zv         = gradNLL dist mSErr xss ys tree theta_t VS.! ix
        inv_slope' = min 4.0 . max 0.0625 . abs $ (tau / (stdErr_i * zv))
        theta_t    = fst
                   $ minimizeNLLWithFixedParam dist mSErr 10 xss ys tree ix 
                   $ theta VS.// [(ix, theta_i + delta * t)]
        nll_cond   = nll dist mSErr xss ys tree theta_t
        tau        = signum delta_t * sqrt (2*nll_cond - 2*nll_opt)
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

predictionCI Profile cov xss tree theta alpha = error "Prediction interval not implemented for profile-t"

createSplines profs (VS.toList -> stdErr) = map spline (zip3 [0..] profs stdErr)
  where
    spline (ix, ProfileT taus thetas _, se)
      | VS.length taus < 2 = (evaluate CSpline [(0, -se), (1, se)], evaluate CSpline [(-se, 0), (se, 1)])
      | otherwise          = (tau2theta, theta2tau)
      where
        tau2theta = evaluate CSpline $ sortOnFirst taus (getRow ix thetas)
        theta2tau = evaluate CSpline $ sortOnFirst (getRow ix thetas) taus
    sortOnFirst xs ys = sortOn fst $ zip (VS.toList xs) (VS.toList ys)
    getRow ix mtx = LA.flatten $ mtx LA.? [ix]

{-
    def splines_sketches(self, tau_scale, p_idx, q_idx):
        '''
        Creates the pairwise splines
        used to calculate the pairwise plot points.

        Returns
        -------
        spline_g : array_like
                    matrix of spline functions for
                    every pair of parameters.
        '''
        spline_g = [[lambda x: x for _ in range(self.likelihood.n)] for _ in range(self.likelihood.n)]
        for p_idx in range(self.likelihood.n):
            for q_idx in range(self.likelihood.n):
                if p_idx == q_idx:
                    continue
                theta_q = self.profiles[q_idx][1][p_idx, :]
                tau_q = self.profiles[q_idx][0]

                gpq = self.spline_theta2tau[p_idx](theta_q)/tau_scale
                idx = np.abs(gpq) < 1
                gpq = np.arccos(gpq[idx])

                if len(tau_q) < 2:
                    spline_g[p_idx][q_idx] = lambda x: x
                else:
                    spline_g[p_idx][q_idx] = CubicSpline(tau_q[idx], gpq)
        return spline_g

    def approximate_contour(self, ix1, ix2, alpha):
        '''
        Approximates de profile countour plot
        for parameters ix1 and ix2 with confidence alpha.

        Parameters
        ----------
        ix1 : array_like
               index of the first parameter
        ix2 : array_like
               index of the second parameter
        alpha : float
                 significance level

        Returns
        -------
        p : array_like
             points for the first parameter
        q : array_like
             points for the second parameter
        '''
        tau_scale = np.sqrt(self.likelihood.n *
                            stats.f.ppf(1 - alpha,
                                        self.likelihood.n,
                                        self.likelihood.m - self.likelihood.n)
                            )
        spline_g = self.splines_sketches(tau_scale)

        angle_pairs = [(0, spline_g[ix2][ix1](1)),
                       (spline_g[ix1][ix2](1), 0),
                       (np.pi, spline_g[ix2][ix1](-1)),
                       (spline_g[ix1][ix2](-1), np.pi)
                       ]

        a = np.zeros(5)
        d = np.zeros(5)
        for j in range(4):
            a_j = (angle_pairs[j][0] + angle_pairs[j][1]) / 2.0
            d_j = angle_pairs[j][0] - angle_pairs[j][1]
            if d_j < 0:
                d_j = -d_j
                a_j = -a_j
            a[j] = a_j
            d[j] = d_j
        a[4] = a[0] + 2*np.pi
        d[4] = d[0]

        ixs = np.argsort(a)
        a = a[ixs]
        d = d[ixs]

        spline_ad = CubicSpline(a, d)
        n_steps = 100
        taup = np.zeros(n_steps)
        tauq = np.zeros(n_steps)
        p = np.zeros(n_steps)
        q = np.zeros(n_steps)
        for i in range(n_steps):
            a_i = i * np.pi * 2 / (n_steps - 1) - np.pi
            d_i = spline_ad(a_i)
            taup[i] = np.cos(a_i + d_i / 2) * tau_scale
            tauq[i] = np.cos(a_i - d_i / 2) * tau_scale
            p[i] = self.spline_tau2theta[ix1](taup[i])
            q[i] = self.spline_tau2theta[ix2](tauq[i])
        return p, q
-}