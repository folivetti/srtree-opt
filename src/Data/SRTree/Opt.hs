{-# language FlexibleInstances #-}
{-# language OverloadedStrings #-}
{-# language ImportQualifiedPost #-}
{-# language ViewPatterns #-}
module Data.SRTree.Opt
    ( optimize, sse, mse, rmse, Column, Columns, minimizeGaussian, minimizeBinomial, nll, Distribution (..) )
    where

import Data.SRTree (SRTree (..), gradParams, floatConstsToParam, gradParams, evalTree)
import Data.SRTree.Recursion (Fix(..))
import Data.Vector qualified as V
import Numeric.GSL.Fitting (FittingMethod (..), nlFitting)
import Numeric.GSL.Minimization (MinimizeMethodD (..), minimizeVD)
import Numeric.LinearAlgebra qualified as LA
import Data.SRTree.Likelihoods

optimize :: Maybe Distribution -> Maybe Double -> Int -> Columns -> Column -> Fix SRTree -> (Fix SRTree, V.Vector Double)
optimize mdist msErr niter xss ys tree = (tree', opt)
  where
    (tree', V.fromList -> t0) = floatConstsToParam tree
    opt         = case mdist of
                    Nothing -> leastSquares niter xss ys tree' t0
                    Just d  -> minimizeNLL d msErr niter xss ys tree' t0

leastSquares :: Int -> Columns -> Column -> Fix SRTree -> V.Vector Double -> V.Vector Double
leastSquares niter xss ys tree t0
  | n == 0    = t0
  | n > m     = t0
  | otherwise = V.fromList . LA.toList . fst $ nlFitting LevenbergMarquardtScaled 1e-6 1e-6 niter model jacob t0'
  where
    t0'         = LA.fromList $ V.toList t0
    n           = LA.size t0'
    m           = LA.size ys
    model theta = let theta' = V.fromList (LA.toList theta)
                   in subtract ys $ evalTree xss theta' LA.scalar tree
    jacob theta = let theta' = V.fromList (LA.toList theta) 
                   in LA.fromColumns . snd $ gradParams xss theta' LA.scalar tree

minimizeNLL :: Distribution -> Maybe Double -> Int -> Columns -> Column -> Fix SRTree -> V.Vector Double -> V.Vector Double
minimizeNLL dist msErr niter xss ys tree t0
  | n == 0    = t0
  | n > m     = t0
  | otherwise = V.fromList . LA.toList . fst $ minimizeVD VectorBFGS2 1e-6 niter 1e-2 1e-6 model jacob t0'
  where
    t0'   = LA.fromList $ V.toList t0
    n     = LA.size t0'
    m     = LA.size ys
    model = nll dist msErr xss ys tree . V.fromList . LA.toList
    jacob = gradNLL dist msErr xss ys tree . V.fromList . LA.toList

minimizeGaussian :: Int -> Columns -> Column -> Fix SRTree -> V.Vector Double -> V.Vector Double
minimizeGaussian = minimizeNLL Gaussian Nothing

minimizeBinomial :: Int -> Columns -> Column -> Fix SRTree -> V.Vector Double -> V.Vector Double
minimizeBinomial = minimizeNLL Bernoulli Nothing
