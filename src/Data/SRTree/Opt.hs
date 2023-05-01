{-# language FlexibleInstances #-}
{-# language OverloadedStrings #-}
{-# language ImportQualifiedPost #-}
module Data.SRTree.Opt
    ( optimize, sse, mse, rmse, Column, Columns )
    where

import Control.Arrow ((&&&))
import Control.Monad.Reader (runReader)
import Control.Monad.State (MonadState (get), State, evalState, gets, modify)
import Data.ByteString.Char8 qualified as B
import Data.ByteString.Lazy qualified as BS
import Data.List (sortOn, isInfixOf)
import Data.Maybe (fromJust)
import Data.SRTree (SRTree (..), gradParams, floatConstsToParam, gradParams, evalTree)
import Data.SRTree.Recursion (Fix(..), cata)
import Data.Vector qualified as V
import Numeric.GSL.Fitting (FittingMethod (..), nlFitting)
import Numeric.LinearAlgebra qualified as LA
import Codec.Compression.GZip
type Columns = V.Vector Column
type Column  = LA.Vector Double

optimize :: Int -> Columns -> Column -> Fix SRTree -> (Fix SRTree, V.Vector Double)
optimize niter xss ys tree = let (tree', t0) = floatConstsToParam tree
                              in (tree', leastSquares niter tree' xss ys $ V.fromList t0)

leastSquares :: Int -> Fix SRTree -> Columns -> Column -> V.Vector Double -> V.Vector Double
leastSquares niter tree xss ys t0
  | n == 0    = t0
  | n > m     = t0
  | otherwise = V.fromList . LA.toList . fst $ nlFitting LevenbergMarquardtScaled 1e-6 1e-6 niter model jacob t0'
  where
    t0'                = LA.fromList $ V.toList t0
    n                  = LA.size t0'
    m                  = LA.size ys
    model theta = let theta' = V.fromList (LA.toList theta) 
                   in subtract ys $ evalTree xss theta' LA.scalar tree
    jacob theta = let theta' = V.fromList (LA.toList theta) 
                   in LA.fromColumns . snd $ gradParams xss theta' LA.scalar tree

sse :: Columns -> Column -> V.Vector Double -> Fix SRTree -> Double
sse xss ys theta tree = sum err
  where
    yhat  = evalTree xss theta LA.scalar tree
    err   = LA.toList $ (ys - yhat) ^ (2 :: Int)

mse :: Columns -> Column -> V.Vector Double -> Fix SRTree -> Double
mse xss ys theta tree = sse xss ys theta tree / fromIntegral (LA.size ys)

rmse :: Columns -> Column -> V.Vector Double -> Fix SRTree -> Double
rmse xss ys theta = sqrt . mse xss ys theta
