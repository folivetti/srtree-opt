{-# language FlexibleInstances #-}
{-# language OverloadedStrings #-}
{-# language ImportQualifiedPost #-}
module Data.SRTree.Opt
    ( loadDataset, optimize, splitTrainVal, sse, mse, rmse, Column, Columns, constToParam, varToConst, paramToVar, getTheta )
    where

import Control.Arrow ((&&&))
import Control.Monad.Reader (runReader)
import Control.Monad.State (MonadState (get), State, evalState, gets, modify)
import Data.ByteString.Char8 qualified as B
import Data.ByteString.Lazy qualified as BS
import Data.List (sortOn, isInfixOf)
import Data.Maybe (fromJust)
import Data.SRTree (OptIntPow (..), SRTree (..), deriveBy, evalTree, evalTreeMap)
import Data.Vector qualified as V
import Numeric.GSL.Fitting (FittingMethod (..), nlFitting)
import Numeric.LinearAlgebra qualified as LA
import Codec.Compression.GZip
type Columns = V.Vector Column
type Column  = LA.Vector Double

byteshow :: Int -> B.ByteString
byteshow = B.pack . show
{-# INLINE byteshow #-}

loadMtx :: [[B.ByteString]] -> LA.Matrix Double
loadMtx = LA.fromLists . map (map (read . B.unpack))
{-# INLINE loadMtx #-}

toColumns :: [Int] -> LA.Matrix Double -> [LA.Vector Double]
toColumns ixs = LA.toColumns . (LA.¿ ixs)
{-# INLINE toColumns #-}

isGZip :: String -> Bool
isGZip [ ]    = False
isGZip [_]    = False
isGZip "gz"   = True
isGZip [_, _] = False
isGZip (_:xs) = isGZip xs
{-# INLINE isGZip #-}

loadDataset :: String -> [Int] -> Int -> Bool -> IO ((Columns, Column), [(B.ByteString, Int)])
loadDataset filename columns target hasHeader = do
  let 
    appGz = if isGZip filename then decompress else id
    separator 
      | "csv" `isInfixOf` filename = ','
      | "tsv" `isInfixOf` filename = '\t'
      | otherwise                  = error "unsupported format"
  content <- filter (not.null) . map (B.split separator) . B.split '\n' . B.pack . map (toEnum . fromEnum) . BS.unpack . appGz <$> BS.readFile filename
  let
    datum     =  loadMtx $ (if hasHeader then tail else id) content
    (_, n)    = LA.size datum
    columns'  = if null columns then [0 .. n - 2] else columns
    target'   = if target < 0 then [n - 1] else [target]
    (xss, ys) = (V.fromList . toColumns columns' &&& head . toColumns target' ) datum
    header    = if hasHeader
                  then sortOn snd $ zip (head content) columns'
                  else map (("x" <>) . byteshow &&& id) [0 .. length columns' - 1]
  pure ((xss, ys), header)

splitTrainVal :: Int -> (Columns, Column) -> ((Columns, Column), (Columns, Column))
splitTrainVal nrows (xss, ys)
  | nrows >= n = error "selecting more rows than the size of the data set"
  | nrows <= 0 = ( (xss, ys), (xss, ys) )
  | otherwise  = ( (V.map getTrain xss, getTrain ys), (V.map getVal xss, getVal ys) )
  where
    n        = LA.size ys
    getTrain = LA.subVector 0 nrows
    getVal   = LA.subVector nrows (n - nrows)

getTheta :: SRTree Int Double -> Column
getTheta = LA.fromList . go
  where
    go :: SRTree Int Double -> [Double]
    go Empty         = []
    go (Var _)       = []
    go (Param _)     = []
    go (Const v)     = [v | fromIntegral (round v :: Int) /= v]
    go (Fun _ t)     = go t
    go (Pow t _)     = go t
    go (Add l r)     = go l <> go r
    go (Sub l r)     = go l <> go r
    go (Mul l r)     = go l <> go r
    go (Div l r)     = go l <> go r
    go (Power l r)   = go l <> go r
    go (LogBase l r) = go l <> go r

replaceTheta :: SRTree Int Double -> Column -> SRTree Int Double
replaceTheta tree theta = go tree `evalState` LA.toList theta
  where
    go :: SRTree Int Double -> State [Double] (SRTree Int Double)
    go Empty         = pure Empty
    go (Var ix)      = pure $ Var ix
    go (Param _)     = do { v <- gets head; modify tail; pure $ Const v }
    go (Const v)     = pure $ Const v
    go (Fun f t)     = Fun f <$> go t
    go (Pow t i)     = (`Pow` i) <$> go t
    go (Add l r)     = do { l' <- go l; r' <- go r; pure $ Add l' r' }
    go (Sub l r)     = do { l' <- go l; r' <- go r; pure $ Sub l' r' }
    go (Mul l r)     = do { l' <- go l; r' <- go r; pure $ Mul l' r' }
    go (Div l r)     = do { l' <- go l; r' <- go r; pure $ Div l' r' }
    go (Power l r)   = do { l' <- go l; r' <- go r; pure $ Power l' r' }
    go (LogBase l r) = do { l' <- go l; r' <- go r; pure $ LogBase l' r' }


constToParam :: SRTree Int Double -> SRTree Int Double
constToParam tree = go tree `evalState` 0
  where
    go :: SRTree Int Double -> State Int (SRTree Int Double)
    go Empty         = pure Empty
    go (Var ix)      = pure $ Var ix
    go (Param ix)    = pure $ Param ix
    go (Const v)     = if fromIntegral (round v :: Int) == v
                         then pure (Const v)
                         else do { ix <- get; modify (+1); pure (Param ix) }
    go (Fun f t)     = Fun f <$> go t
    go (Pow t i)     = (`Pow` i) <$> go t
    go (Add l r)     = do { l' <- go l; r' <- go r; pure $ Add l' r' }
    go (Sub l r)     = do { l' <- go l; r' <- go r; pure $ Sub l' r' }
    go (Mul l r)     = do { l' <- go l; r' <- go r; pure $ Mul l' r' }
    go (Div l r)     = do { l' <- go l; r' <- go r; pure $ Div l' r' }
    go (Power l r)   = do { l' <- go l; r' <- go r; pure $ Power l' r' }
    go (LogBase l r) = do { l' <- go l; r' <- go r; pure $ LogBase l' r' }

paramToVar :: SRTree Int a -> SRTree Int a
paramToVar = go
  where
    go :: SRTree Int a -> SRTree Int a
    go Empty         = Empty
    go (Var ix)      = Var ix
    go (Param ix)    = Var ix
    go (Const v)     = Const v
    go (Fun f t)     = Fun f $ go t
    go (Pow t i)     = (`Pow` i) $ go t
    go (Add l r)     = Add (go l) (go r)
    go (Sub l r)     = Sub (go l) (go r)
    go (Mul l r)     = Mul (go l) (go r)
    go (Div l r)     = Div (go l) (go r)
    go (Power l r)   = Power (go l) (go r)
    go (LogBase l r) = LogBase (go l) (go r)

varToConst :: Columns -> SRTree Int Double -> SRTree Int Column
varToConst xss = go
  where
    go :: SRTree Int Double -> SRTree Int Column
    go Empty         = Empty
    go (Var ix)      = Const (xss V.! ix)
    go (Param ix)    = Param ix
    go (Const v)     = Const $ LA.scalar v
    go (Fun f t)     = Fun f $ go t
    go (Pow t i)     = (`Pow` i) $ go t
    go (Add l r)     = Add (go l) (go r)
    go (Sub l r)     = Sub (go l) (go r)
    go (Mul l r)     = Mul (go l) (go r)
    go (Div l r)     = Div (go l) (go r)
    go (Power l r)   = Power (go l) (go r)
    go (LogBase l r) = LogBase (go l) (go r)

optimize :: Int -> Columns -> Column -> SRTree Int Double -> SRTree Int Double
optimize niter xss ys tree = let t0    = getTheta tree
                                 tree' = constToParam tree
                              in replaceTheta tree' $ leastSquares niter tree' xss ys t0

leastSquares :: Int -> SRTree Int Double -> Columns -> Column -> Column -> Column
leastSquares niter tree xss ys t0
  | LA.size t0 == 0 = t0
  | otherwise       = fst $ nlFitting LevenbergMarquardt 1e-6 1e-6 niter model jacob t0
  where
    n                  = LA.size t0
    m                  = LA.size ys
    tree'              = paramToVar $ varToConst xss tree
    derive             = flip deriveBy
    evalSRTree t theta = fromJust $ evalTree t `runReader` (Just . LA.scalar . (theta LA.!))
    model              = subtract ys . evalSRTree tree'

    jacob t | sz == (1,1) = LA.asColumn $ LA.konst (head $ LA.toList $ head j) m
            | otherwise   = LA.fromColumns j
            where
                j  = map ((`evalSRTree` t) . derive tree') [0 .. n-1]
                sz = LA.size $ LA.fromColumns j

sse :: Columns -> Column -> SRTree Int Double -> Double
sse xss ys tree = sum err
  where
    yhat  = fromJust $ evalTreeMap LA.scalar tree `runReader` (xss V.!?)
    err   = LA.toList $ (ys - yhat) ^ (2 :: Int)

mse :: Columns -> Column -> SRTree Int Double -> Double
mse xss ys tree = sse xss ys tree / fromIntegral (LA.size ys)

rmse :: Columns -> Column -> SRTree Int Double -> Double
rmse xss ys = sqrt . mse xss ys

instance OptIntPow (LA.Vector Double) where
  (^.) = (^^)
