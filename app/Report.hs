module Report where

import qualified Data.Vector.Storable as VS
import qualified Numeric.LinearAlgebra as LA
import Data.Maybe ( fromMaybe )
import Statistics.Distribution.FDistribution ( fDistribution )
import Statistics.Distribution ( quantile )
import System.Random ( StdGen, split )
import Data.Random.Normal ( normals )

import Data.SRTree ( SRTree, Fix (..), floatConstsToParam, paramsToConst, countNodes )
import Data.SRTree.AD ( reverseModeUnique )
import Data.SRTree.Likelihoods
import Data.SRTree.ModelSelection ( aic, bic, evidence, logFunctional, logParameters, mdl, mdlFreq, mdlLatt )
import Data.SRTree.ConfidenceIntervals
import Data.SRTree.Opt (optimize, minimizeNLLWithFixedParam, minimizeNLL)
import Data.SRTree.Datasets ( loadDataset )
import Debug.Trace ( trace )

import Args

-- store the datasets split into training, validation and test
data Datasets = DS { _xTr  :: Columns
                   , _yTr  :: Column
                   , _xVal :: Maybe Columns
                   , _yVal :: Maybe Column
                   , _xTe  :: Maybe Columns
                   , _yTe  :: Maybe Column
                   }

-- basic fields name
basicFields :: [String]
basicFields = [ "Index"
              , "Filename"
              , "Expression"
              , "Number_of_nodes"
              , "Number_of_parameters"
              , "Parameters"
              , "Iterations"
              ]

-- basic information about the tree
data BasicInfo = Basic { _index   :: Int
                       , _fname   :: String
                       , _expr    :: Fix SRTree
                       , _nNodes  :: Int
                       , _nParams :: Int
                       , _params  :: [Double]
                       , _iters   :: Int
                       }

-- optimization fields
optFields :: [String]
optFields = [ "SSE_train_orig"
            , "SSE_val_orig"
            , "SSE_test_orig"
            , "SSE_train_opt"
            , "SSE_val_opt"
            , "SSE_test_opt"
            ]

-- optimization information
data SSE = SSE { _sseTr  :: Double
               , _sseVal :: Double
               , _sseTe  :: Double
               }

-- model selection fields
modelFields :: [String]
modelFields = [ "BIC"
              , "BIC_val"
              , "AIC"
              , "AIC_val"
              , "Evidence"
              , "EvidenceVal"
              , "MDL"
              , "MDL_Freq"
              , "MDL_Lattice"
              , "MDL_val"
              , "MDL_Freq_val"
              , "MDL_Lattice_val"
              , "NegLogLikelihood_train"
              , "NegLogLikelihood_val"
              , "NegLogLikelihood_test"
              , "LogFunctional"
              , "LogParameters"
              , "Fisher"
              ]

-- model selection information
data Info = Info { _bic     :: Double
                 , _bicVal  :: Double
                 , _aic     :: Double
                 , _aicVal  :: Double
                 , _evidence :: Double
                 , _evidenceVal :: Double
                 , _mdl     :: Double
                 , _mdlFreq :: Double
                 , _mdlLatt :: Double
                 , _mdlVal  :: Double
                 , _mdlFreqVal :: Double
                 , _mdlLattVal :: Double
                 , _nllTr   :: Double
                 , _nllVal  :: Double
                 , _nllTe   :: Double
                 , _cc      :: Double
                 , _cp      :: Double
                 , _fisher  :: [Double]
                 }

-- load the datasets
getDataset :: Args -> IO (Datasets, String)
getDataset args = do
  ((xTr, yTr, xVal, yVal), varnames) <- loadDataset (dataset args) (hasHeader args)
  let (mXVal, mYVal) = if yTr == yVal
                         then (Nothing, Nothing)
                         else (Just xVal, Just yVal)
  (mXTe, mYTe) <- if null (test args)
                    then pure (Nothing, Nothing)
                    else do ((xTe, yTe, _, _), _) <- loadDataset (test args) (hasHeader args)
                            pure (Just xTe, Just yTe)
  pure (DS xTr yTr mXVal mYVal mXTe mYTe, varnames)

getBasicStats :: Args -> StdGen -> Datasets -> Fix SRTree -> Int -> BasicInfo
getBasicStats args seed dset tree ix
  | anyNaN    = getBasicStats args (snd $ split seed) dset tree ix
  | otherwise = Basic ix (infile args) tOpt nNodes nParams params n
  where
    (tree', theta0) = floatConstsToParam tree
    thetas          = if restart args
                        then take nParams (normals seed)
                        else theta0
    (_, t, n)       = optimize (dist args) (msErr args) (niter args) (_xTr dset) (_yTr dset) (Just thetas) tree
    tOpt            = paramsToConst (VS.toList t) tree'
    nNodes          = countNodes tOpt :: Int
    nParams         = length theta0
    params          = VS.toList t
    anyNaN          = VS.any isNaN t

getSSE :: Datasets -> Fix SRTree -> SSE
getSSE dset tree = SSE tr val te
  where
    (t, th) = floatConstsToParam tree
    tr  = sse (_xTr dset) (_yTr dset) t (VS.fromList th)
    val = case (_xVal dset, _yVal dset) of
            (Nothing, _)           -> 0.0
            (_, Nothing)           -> 0.0
            (Just xVal, Just yVal) -> sse xVal yVal t (VS.fromList th)
    te  = case (_xTe dset, _yTe dset) of
            (Nothing, _)           -> 0.0
            (_, Nothing)           -> 0.0
            (Just xTe, Just yTe)   -> sse xTe yTe t (VS.fromList th)

getInfo :: Args -> Datasets -> Fix SRTree -> Fix SRTree -> Info
getInfo args dset tree treeVal =
  Info { _bic     = bic dist' msErr' xTr yTr thetaOpt' tOpt
       , _bicVal  = bicVal
       , _aic     = aic dist' msErr' xTr yTr thetaOpt' tOpt
       , _aicVal  = aicVal
       , _evidence = evidence dist' msErr' xTr yTr thetaOpt' tOpt
       , _evidenceVal = evidenceVal
       , _mdl     = mdl dist' msErr' xTr yTr thetaOpt' tOpt
       , _mdlFreq = mdlFreq dist' msErr' xTr yTr thetaOpt' tOpt
       , _mdlLatt = mdlLatt dist' msErr' xTr yTr thetaOpt' tOpt
       , _mdlVal  = mdlVal
       , _mdlFreqVal = mdlFreqVal
       , _mdlLattVal = mdlLattVal
       , _nllTr   = nllTr
       , _nllVal  = nllVal
       , _nllTe   = nllTe
       , _cc      = logFunctional tOpt
       , _cp      = logParameters dist' msErr' xTr yTr thetaOpt' tOpt
       , _fisher  = fisherNLL dist' (msErr args) xTr yTr tOpt thetaOpt'
       }
  where
    (xTr, yTr)       = (_xTr dset, _yTr dset)
    (xVal, yVal)     = case (_xVal dset, _yVal dset) of
                         (Nothing, _)     -> (xTr, yTr)
                         (_, Nothing)     -> (xTr, yTr)
                         (Just a, Just b) -> (a, b)
    (tOpt, thetaOpt) = floatConstsToParam tree
    thetaOpt'        = VS.fromList thetaOpt

    (tOptVal, thetaOptVal) = floatConstsToParam treeVal
    thetaOptVal'           = VS.fromList thetaOptVal

    dist'            = fromMaybe Gaussian (dist args)
    msErr'           = msErr args
    nllTr            = nll dist' msErr' (_xTr dset) (_yTr dset) tOpt (VS.fromList thetaOpt)
    bicVal           = case (_xVal dset, _yVal dset) of
                         (Nothing, _) -> 0.0
                         (_, Nothing) -> 0.0
                         _            -> bic dist' msErr' xVal yVal thetaOptVal' tOptVal
    aicVal           = case (_xVal dset, _yVal dset) of
                         (Nothing, _) -> 0.0
                         (_, Nothing) -> 0.0
                         _            -> aic dist' msErr' xVal yVal thetaOptVal' tOptVal
    evidenceVal      = case (_xVal dset, _yVal dset) of
                         (Nothing, _) -> 0.0
                         (_, Nothing) -> 0.0
                         _            -> evidence dist' msErr' xVal yVal thetaOptVal' tOptVal
    nllVal           = case (_xVal dset, _yVal dset) of
                         (Nothing, _) -> 0.0
                         (_, Nothing) -> 0.0
                         _            -> nll dist' msErr' xVal yVal tOptVal (VS.fromList thetaOptVal)
    mdlVal           = case (_xVal dset, _yVal dset) of
                         (Nothing, _) -> 0.0
                         (_, Nothing) -> 0.0
                         _            -> mdl dist' msErr' xVal yVal thetaOptVal' tOptVal
    mdlFreqVal       = case (_xVal dset, _yVal dset) of
                         (Nothing, _) -> 0.0
                         (_, Nothing) -> 0.0
                         _            -> mdlFreq dist' msErr' xVal yVal thetaOptVal' tOptVal
    mdlLattVal       = case (_xVal dset, _yVal dset) of
                         (Nothing, _) -> 0.0
                         (_, Nothing) -> 0.0
                         _            -> mdlLatt dist' msErr' xVal yVal thetaOptVal' tOptVal
    nllTe            = case (_xTe dset, _yTe dset) of
                         (Nothing, _)           -> 0.0
                         (_, Nothing)           -> 0.0
                         (Just xTe, Just yTe) -> nll dist' msErr' xTe yTe tOpt (VS.fromList thetaOpt)

getCI :: Args -> Datasets -> BasicInfo -> Double -> (BasicStats, [CI], [CI], [CI], [CI])
getCI args dset basic alpha' = (stats', cis, pis_tr, pis_val, pis_te)
  where
    (tree, _)  = floatConstsToParam (_expr basic)
    theta      = _params basic
    tau_max    = (quantile (fDistribution (_nParams basic) (LA.size yTr - _nParams basic)) (1 - 0.01))
    (xTr, yTr) = (_xTr dset, _yTr dset)
    dist'      = fromMaybe Gaussian (dist args)
    msErr'     = msErr args
    stats'     = getStatsFromModel dist' msErr' xTr yTr tree (VS.fromList theta)
    profiles   = getAllProfiles dist' msErr' xTr yTr tree (VS.fromList theta) (_stdErr stats')
    method     = if useProfile args
                   then Profile stats' profiles
                   else Laplace stats'
    predFun   = predict dist' tree (VS.fromList theta)
    stdErr    = _stdErr stats' VS.! 0

    prof th t = let (thOpt, _) = minimizeNLL dist' msErr' 1000 xTr yTr t th
                 in case getProfile dist' msErr' xTr yTr t thOpt stdErr tau_max 0 of
                      Left th' -> prof th' t
                      Right p  -> (_tau2theta p, _opt p)
    jac xss   = LA.toRows . LA.fromColumns $ snd $ reverseModeUnique xss (VS.fromList theta) VS.singleton tree

    cis       = paramCI method (LA.size yTr) (VS.fromList theta) alpha' 
    pis_tr    = predictionCI method dist' predFun jac prof xTr tree (VS.fromList theta) alpha' 
    pis_val   = case (_xVal dset, _yVal dset) of
                  (Nothing, _)   -> []
                  (Just xVal, _) -> predictionCI method dist' predFun jac prof xVal tree (VS.fromList theta) alpha'
    pis_te    = case (_xTe dset, _yTe dset) of
                  (Nothing, _)  -> []
                  (Just xTe, _) -> predictionCI method dist' predFun jac prof xTe tree (VS.fromList theta) alpha'
