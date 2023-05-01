{-# language LambdaCase #-}
{-# language ImportQualifiedPost #-}
module Main (main) where

import Control.Monad ( unless, forM_ )
import Data.Char (toLower, toUpper)
import Data.List (intercalate)
import Data.SRTree ( SRTree )
import Data.SRTree.Recursion ( Fix )
import Data.SRTree.Opt (optimize, sse)
import Data.SRTree.Print qualified as P
import Data.SRTree.Datasets (loadDataset)
import Options.Applicative
import System.IO ( hClose, hPutStrLn, openFile, stderr, stdout, IOMode(WriteMode), Handle )
import Text.ParseSR (SRAlgs (..))
import Text.ParseSR.IO (withInput)
import Text.Read (readMaybe)

envelope :: a -> [a] -> [a]
envelope c xs = c : xs <> [c]
{-# INLINE envelope #-}

sralgsHelp :: [String]
sralgsHelp = map (envelope '\'' . map toLower . show) [toEnum 0 :: SRAlgs ..]
{-# INLINE sralgsHelp #-}

sralgsReader :: ReadM SRAlgs
sralgsReader = do
  sr <- map toUpper <$> str
  eitherReader $ case readMaybe sr of
    Nothing -> pure . Left $ "unknown algorithm. Available options are " <> intercalate "," sralgsHelp
    Just x  -> pure . Right $ x

columnsReader :: ReadM [Int]
columnsReader = do
  colsStr <- ('[':) . (<> "]") <$> str
  eitherReader $ case readMaybe colsStr of
    Nothing -> pure . Left $ "wrong format " <> colsStr
    Just x  -> pure . Right $ x

data Args = Args
    {   from        :: SRAlgs
      , infile      :: String
      , outfile     :: String
      , stats       :: String
      , dataset     :: String
      , niter       :: Int
      , hasHeader   :: Bool
      , simpl       :: Bool
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

openWriteWithDefault :: Handle -> String -> IO Handle
openWriteWithDefault dflt fname = 
    if null fname 
       then pure dflt 
       else openFile fname WriteMode
{-# INLINE openWriteWithDefault #-}

printResults :: String -> String -> (Fix SRTree -> (Fix SRTree, String)) -> [Either String (Fix SRTree)] -> IO ()
printResults fname sname f exprs = do
  hExpr <- openWriteWithDefault stdout fname
  hStat <- openWriteWithDefault stderr sname
  hPutStrLn hStat "SSE_train_orig,SSE_val_orig,SSE_train_opt,SSE_val_opt"
  forM_ exprs $ \case 
                   Left  err -> do hPutStrLn hExpr $ "invalid expression: " <> err
                                   hPutStrLn hStat $ "invalid expression: " <> err
                   Right ex  -> do let (ex', sts) = f ex
                                   hPutStrLn hExpr (P.showExpr ex')
                                   hPutStrLn hStat sts
  unless (null fname) $ hClose hExpr
  unless (null sname) $ hClose hStat

main :: IO ()
main = do
  args <- execParser opts
  ((xTr, yTr, xVal, yVal), varnames) <- loadDataset (dataset args) (hasHeader args)
  let optimizer = optimize (niter args) xTr yTr
      sseTr t    = show . sse xTr yTr t
      sseVal t   = show . sse xVal yVal t
      genStats  tree = let (t, theta) = optimizer tree
                        in (t, intercalate "," [sseTr theta tree, sseVal theta tree, sseTr theta t, sseVal theta t])
  withInput (infile args) (from args) varnames False (simpl args)
    >>= printResults (outfile args) (stats args) genStats
  
  where 
      opts = info (opt <**> helper)
            ( fullDesc <> progDesc "Optimize the parameters of Symbolic Regression expressions."
            <> header "srtree-opt - a CLI tool to (re)optimize the numeric parameters of symbolic regression expressions"
            )
