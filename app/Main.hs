{-# language LambdaCase #-}
{-# language ImportQualifiedPost #-}
module Main (main) where

import Control.Monad ( unless, forM_ )
import Data.Bifunctor (first)
import Data.ByteString.Char8 qualified as B
import Data.Char (toLower, toUpper)
import Data.List (intercalate)
import Data.SRTree ( SRTree )
import Data.SRTree.EqSat (simplifyEqSat)
import Data.SRTree.Opt (Column, Columns, loadDataset, optimize, splitTrainVal, sse)
import Data.SRTree.Print qualified as P
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
      , trainRows   :: Int
      , cols        :: [Int]
      , target      :: Int
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
       <> help "Filename of the dataset used for optimizing the parameters." )
   <*> option auto
       ( long "rows"
       <> short 'r'
       <> metavar "ROWS"
       <> showDefault
       <> value 0
       <> help "Number of rows to use as training data, the remainder will be used as validation. Values <= 0 will use the whole data as training and validation.")
   <*> option columnsReader
       ( long "columns"
       <> short 'c'
       <> metavar "COLUMNS"
       <> showDefault
       <> value []
       <> help "Index of columns to use as variables. Default \"\" uses all but the last column.")
   <*> option auto
       ( long "target"
       <> short 't'
       <> metavar "TARGET"
       <> showDefault
       <> value (-1)
       <> help "Index of colum to use as the target variable. Default (-1) uses the last column.")
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

openData :: Args -> IO (((Columns, Column), (Columns, Column)), [(B.ByteString, Int)])
openData args = first (splitTrainVal (trainRows args)) 
             <$> loadDataset (dataset args  ) (cols args) (target args) (hasHeader args)

openWriteWithDefault :: Handle -> String -> IO Handle
openWriteWithDefault dflt fname = 
    if null fname 
       then pure dflt 
       else openFile fname WriteMode
{-# INLINE openWriteWithDefault #-}

printResults :: String -> String -> (SRTree Int Double -> (SRTree Int Double, String)) -> [Either String (SRTree Int Double)] -> IO ()
printResults fname sname f exprs = do
  hExpr <- openWriteWithDefault stdout fname
  hStat <- openWriteWithDefault stderr sname
  hPutStrLn hStat "SSE_train_orig,SSE_val_orig,SSE_train_opt,SSE_val_opt"
  forM_ exprs $ \case 
                   Left  err -> do hPutStrLn hExpr $ "invalid expression: " <> err
                                   hPutStrLn hStat $ "invalid expression: " <> err
                   Right ex  -> do let (ex', sts) = f ex
                                   hPutStrLn hExpr (P.showDefault ex')
                                   hPutStrLn hStat sts
  unless (null fname) $ hClose hExpr
  unless (null sname) $ hClose hStat

main :: IO ()
main = do
  args <- execParser opts
  (((xTr, yTr),(xVal, yVal)), headers) <- openData args
  let optimizer = optimize (niter args) xTr yTr
      varnames  = intercalate "," (map (B.unpack.fst) headers)
      sseTr     = show . sse xTr yTr
      sseVal    = show . sse xVal yVal
      genStats  tree = let tree' = if simpl args then simplifyEqSat tree else tree
                           t = optimizer tree'
                        in (t, intercalate "," [sseTr tree', sseVal tree', sseTr t, sseVal t])
  withInput (infile args) (from args) varnames False
    >>= printResults (outfile args) (stats args) genStats
  
  where 
      opts = info (opt <**> helper)
            ( fullDesc <> progDesc "Optimize the parameters of Symbolic Regression expressions."
            <> header "srtree-opt - a CLI tool to (re)optimize the numeric parameters of symbolic regression expressions"
            )
