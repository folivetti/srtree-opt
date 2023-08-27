{-# language LambdaCase #-}
{-# language ImportQualifiedPost #-}
{-# language ViewPatterns #-}
module Main (main) where

import Control.Monad ( unless, forM_ )
import Data.Char (toLower, toUpper)
import Data.List (intercalate)
import Data.Maybe ( fromMaybe )
import Data.SRTree ( SRTree, Fix (..), floatConstsToParam, paramsToConst, countNodes )
import Data.SRTree.Opt (optimize, sse, Distribution (..), nll)
import Data.SRTree.AD
import Data.SRTree.Likelihoods
import Data.SRTree.ModelSelection
import Data.SRTree.Print qualified as P
import Data.SRTree.Datasets (loadDataset)
import Data.SRTree.ConfidenceIntervals
import Options.Applicative
import System.IO ( hClose, hPutStrLn, openFile, stderr, stdout, IOMode(WriteMode), Handle )
import Text.ParseSR (SRAlgs (..))
import Text.ParseSR.IO (withInput)
import Text.Read (readMaybe)
import Data.Vector.Storable qualified as VS
import System.Random
import Data.Random.Normal ( normals )
import qualified Numeric.LinearAlgebra as LA 

import Debug.Trace ( trace )

envelope :: a -> [a] -> [a]
envelope c xs = c : xs <> [c]
{-# INLINE envelope #-}

sralgsHelp :: [String]
sralgsHelp = map (envelope '\'' . map toLower . show) [toEnum 0 :: SRAlgs ..]
{-# INLINE sralgsHelp #-}

distHelp :: [String]
distHelp = map (envelope '\'' . map toLower . show) [toEnum 0 :: Distribution ..]
{-# INLINE distHelp #-}

sralgsReader :: ReadM SRAlgs
sralgsReader = do
  sr <- map toUpper <$> str
  eitherReader $ case readMaybe sr of
    Nothing -> pure . Left $ "unknown algorithm. Available options are " <> intercalate "," sralgsHelp
    Just x  -> pure . Right $ x

s2Reader :: ReadM (Maybe Double)
s2Reader = do
  s <- str
  eitherReader $ case readMaybe s of
    Nothing -> pure . Left $ "wrong format " <> s
    mx      -> pure . Right $ mx

distRead :: ReadM (Maybe Distribution)
distRead = do
  d <- capitalize <$> str
  eitherReader $ case readMaybe d of
                    Nothing -> pure . Left $ "unsupported distribution " <> d
                    Just x  -> pure . Right $ Just x
  where
    capitalize ""     = ""
    capitalize (c:cs) = toUpper c : map toLower cs

data Args = Args
    {   from        :: SRAlgs
      , infile      :: String
      , outfile     :: String
      , stats       :: String
      , dataset     :: String
      , test        :: String
      , niter       :: Int
      , hasHeader   :: Bool
      , simpl       :: Bool
      , dist        :: Maybe Distribution
      , msErr       :: Maybe Double
      , restart     :: Bool
      , rseed       :: Int
      , toScreen    :: Bool
      , useProfile  :: Bool
    } deriving Show

opt :: Parser Args
opt = Args
   <$> option sralgsReader
       ( long "from"
       <> short 'f'
       <> metavar ("[" <> intercalate "|" sralgsHelp <> "]")
       <> help "Input expression format" )
   <*> strOption
       ( long "input"
       <> short 'i'
       <> metavar "INPUT-FILE"
       <> showDefault
       <> value ""
       <> help "Input file containing expressions. Empty string gets expression from stdin." )
   <*> strOption
       ( long "output"
       <> short 'o'
       <> metavar "OUTPUT-FILE"
       <> showDefault
       <> value ""
       <> help "Output file to store expressions. Empty string prints expressions to stdout." )
   <*> strOption
       ( long "stats"
       <> short 's'
       <> metavar "STATS-FILE"
       <> showDefault
       <> value ""
       <> help "Output file to store the sse of the trees. Empty string prints expressions to stderr." )
   <*> strOption
       ( long "dataset"
       <> short 'd'
       <> metavar "DATASET-FILENAME"
       <> help "Filename of the dataset used for optimizing the parameters. Empty string omits stats that make use of the training data. It will auto-detect and handle gzipped file based on gz extension. It will also auto-detect the delimiter. \nThe filename can include extra information: filename.csv:start:end:target:vars where start and end corresponds to the range of rows that should be used for fitting, target is the column index (or name) of the target variable and cols is a comma separated list of column indeces or names of the variables in the same order as used by the symbolic model." )
   <*> strOption
       ( long "test"
       <> metavar "TEST"
       <> showDefault
       <> value ""
       <> help "Filename of the test dataset. Empty string omits stats that make use of the training data. It can have additional information as in the training set, but the validation range will be discarded." )
   <*> option auto
       ( long "niter"
       <> metavar "NITER"
       <> showDefault
       <> value 10
       <> help "Number of iterations for the optimization algorithm.")
   <*> switch
       ( long "hasheader"
       <> help "Uses the first row of the csv file as header.")
   <*> switch
        ( long "simplify"
        <> help "Apply basic simplification." )
   <*> option distRead
        ( long "distribution"
        <> metavar ("[" <> intercalate "|" distHelp <> "]")
        <> showDefault
        <> value Nothing
        <> help "Minimize negative log-likelihood following one of the avaliable distributions. The default will use least squares to optimize the model."
        )
   <*> option s2Reader
       ( long "sErr"
       <> metavar "Serr"
       <> showDefault
       <> value Nothing
       <> help "Estimated standard error of the data. If not passed, it uses the model MSE.")
   <*> switch
        ( long "restart"
        <> help "If set, it samples the initial values of the parameters using a Gaussian distribution N(0, 1), otherwise it uses the original values of the expression." )
   <*> option auto
       ( long "seed"
       <> metavar "SEED"
       <> showDefault
       <> value (-1)
       <> help "Random seed to initialize the parameters values. Used only if restart is enabled.")
   <*> switch
        ( long "report"
        <> help "If set, reports the analysis in a user-friendly format instead of csv. It will also include confidence interval for the parameters and predictions" )
   <*> switch
        ( long "profile"
        <> help "If set, it will use profile likelihood to calculate the CIs." )

openWriteWithDefault :: Handle -> String -> IO Handle
openWriteWithDefault dflt fname =
    if null fname
       then pure dflt
       else openFile fname WriteMode
{-# INLINE openWriteWithDefault #-}

csvHeader :: String
csvHeader = "Index,Filename,Expression,Number_of_nodes,Number_of_parameters,Parameters," -- basic
         <> "Iterations,SSE_train_orig,SSE_val_orig,SSE_test_orig,SSE_train_opt,SSE_val_opt,SSE_test_opt," -- sse
         <> "BIC,AIC,MDL,MDL_Freq,NegLogLikelihood_train,NegLogLikelihood_val,NegLogLikelihood_test,LogFunctional,LogParameters,Fisher" -- info
{-# inline csvHeader #-}

printResults :: String -> String -> (Int -> Fix SRTree -> (Fix SRTree, [String])) -> [Either String (Fix SRTree)] -> IO ()
printResults fname sname f exprs = do
  hExpr <- openWriteWithDefault stdout fname
  hStat <- openWriteWithDefault stderr sname
  hPutStrLn hStat csvHeader 
  forM_ (zip [0..] exprs) $ \(ix, e) -> case e of
                   Left  err -> do hPutStrLn hExpr $ "invalid expression: " <> err
                                   hPutStrLn hStat $ "invalid expression: " <> err
                   Right ex  -> do let (ex', sts) = f ix ex
                                   hPutStrLn hExpr (P.showExpr ex')
                                   hPutStrLn hStat (intercalate "," sts)
  unless (null fname) $ hClose hExpr
  unless (null sname) $ hClose hStat

-- printResultsScreen :: String -> String -> (Int -> Fix SRTree -> (Fix SRTree, [String])) -> [Either String (Fix SRTree)] -> IO ()
printResultsScreen f g exprs = do
  forM_ (zip [0..] exprs) $ 
      \(ix, e) -> case e of
                   Left  err -> do putStrLn $ "invalid expression: " <> err
                   Right ex  -> do let (ex', sts) = f ix ex
                                       (sts', cis, pis1, pis2, pis3) = g ex'
                                   -- mapM_ putStrLn sts
                                   print sts'
                                   mapM_ print cis
                                   mapM_ print pis1
                                   mapM_ print pis2
                                   mapM_ print pis3

getInfoFromTree fname fTheta optimizer fisherFun dist' msErr' xTr yTr xVal yVal xTe yTe ix tree = (t', statsStr)
  where
    (tOpt, thetaOpt, its) = optimizer (Just thetas) tree

    t0       = snd $ floatConstsToParam tree
    thetas   = fTheta t0
    t'       = paramsToConst (VS.toList thetaOpt) tOpt
    nNodes   = countNodes t' :: Int
    nParams  = VS.length thetaOpt
    params   = intercalate ";" . map show $ VS.toList thetaOpt
    fisher   = intercalate ";" . map show $ fisherFun tOpt thetaOpt
    statsStr = [ show ix, fname
               , P.showExpr t'
               , show nNodes
               , show nParams
               , params
               , show its
               , show $ sse xTr yTr tree (VS.fromList t0)
               , show $ sse xVal yVal tree (VS.fromList t0)
               , show $ sse xTe yTe tree (VS.fromList t0)
               , show $ sse xTr yTr tOpt thetaOpt
               , show $ sse xVal yVal tOpt thetaOpt
               , show $ sse xTe yTe tOpt thetaOpt
               , show $ bic dist' msErr' xTr yTr thetaOpt tOpt
               , show $ aic dist' msErr' xTr yTr thetaOpt tOpt
               , show $ mdl dist' msErr' xTr yTr thetaOpt tOpt
               , show $ mdlFreq dist' msErr' xTr yTr thetaOpt tOpt
               , show $ nll dist' msErr' xTr yTr tOpt thetaOpt
               , show $ nll dist' msErr' xVal yVal tOpt thetaOpt
               , show $ nll dist' msErr' xTe yTe tOpt thetaOpt
               , show $ logFunctional tOpt
               , show $ logParameters dist' msErr' xTr yTr thetaOpt tOpt
               , fisher
               ]

getStatsFromTree optimizer fisherFun dist' msErr' xTr yTr xVal yVal xTe yTe tree = (stats', cis, pis_tr, pis_val, pis_te)
  where
    (tree', theta) = floatConstsToParam tree
    theta'         = VS.fromList theta
    stats'         = getStatsFromModel dist' msErr' xTr yTr tree' theta'
    predFun        = predict dist' tree' theta'
    jacFun xss     = LA.toRows . LA.fromColumns $ snd $ reverseModeUnique xss theta' VS.singleton tree'
    profiles       = getAllProfiles dist' msErr' xTr yTr tree' theta' (_stdErr stats')
    cis            = paramCI (Profile stats' profiles) (LA.size yTr) theta' 0.05
    pis_tr         = predictionCI (Laplace stats') predFun jacFun (const id) xTr tree' theta' 0.05 
    pis_val        = predictionCI (Laplace stats') predFun jacFun (const id) xVal tree' theta' 0.05 
    pis_te         = predictionCI (Laplace stats') predFun jacFun (const id) xTe tree' theta' 0.05 

--toCSV :: IO ()
toCSV args = do
  g    <- getStdGen
  ((xTr, yTr, xVal, yVal), varnames) <- loadDataset (dataset args) (hasHeader args)
  ((xTe, yTe, _, _), _)              <- if null (test args)
                                          then pure ((xVal, yVal, xVal, yVal), varnames)
                                          else loadDataset (test args) (hasHeader args)
  let seed      = if rseed args < 0 then g else mkStdGen (rseed args)
      optimizer = optimize (dist args) (msErr args) (niter args) xTr yTr
      eTr   t   = show . sse xTr yTr t
      eVal  t   = show . sse xVal yVal t
      eTest t   = show . sse xTe yTe t
      dist'     = fromMaybe Gaussian (dist args)
      fTheta t0 = if restart args
                     then take (length t0) (normals seed)
                     else t0
      fisherFun = fisherNLL dist' (msErr args) xTr yTr
      genStats  = getInfoFromTree (infile args) fTheta optimizer fisherFun dist' (msErr args) xTr yTr xVal yVal xTe yTe

  withInput (infile args) (from args) varnames False (simpl args)
    >>= printResults (outfile args) (stats args) genStats

reportResults args = do
  g    <- getStdGen
  ((xTr, yTr, xVal, yVal), varnames) <- loadDataset (dataset args) (hasHeader args)
  ((xTe, yTe, _, _), _)              <- if null (test args)
                                          then pure ((xVal, yVal, xVal, yVal), varnames)
                                          else loadDataset (test args) (hasHeader args)
  let seed      = if rseed args < 0 then g else mkStdGen (rseed args)
      optimizer = optimize (dist args) (msErr args) (niter args) xTr yTr
      eTr   t   = show . sse xTr yTr t
      eVal  t   = show . sse xVal yVal t
      eTest t   = show . sse xTe yTe t
      dist'     = fromMaybe Gaussian (dist args)
      fTheta t0 = if restart args
                     then take (length t0) (normals seed)
                     else t0
      fisherFun = fisherNLL dist' (msErr args) xTr yTr
      genStats  = getInfoFromTree (infile args) fTheta optimizer fisherFun dist' (msErr args) xTr yTr xVal yVal xTe yTe
      genCIS = getStatsFromTree optimizer fisherFun dist' (msErr args) xTr yTr xVal yVal xTe yTe

  withInput (infile args) (from args) varnames False (simpl args)
    >>= printResultsScreen genStats genCIS

main :: IO ()
main = do
  args <- execParser opts
  if toScreen args
    then reportResults args
    else toCSV args
  where
      opts = info (opt <**> helper)
            ( fullDesc <> progDesc "Optimize the parameters of Symbolic Regression expressions."
            <> header "srtree-opt - a CLI tool to (re)optimize the numeric parameters of symbolic regression expressions"
            )
