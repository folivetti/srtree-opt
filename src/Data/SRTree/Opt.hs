{-# language FlexibleInstances #-}
{-# language OverloadedStrings #-}
{-# language ImportQualifiedPost #-}
{-# language ViewPatterns #-}
module Data.SRTree.Opt
    ( optimize, sse, mse, rmse, Column, Columns, minimizeGaussian, minimizeBinomial, nll, Distribution (..) )
    where

import Data.SRTree (SRTree (..), gradParamsRev, floatConstsToParam, evalTree)
import Data.SRTree.Recursion (Fix(..))
import Data.Vector qualified as V
import Numeric.GSL.Fitting (FittingMethod (..), nlFitting)
import Numeric.GSL.Minimization (MinimizeMethodD (..), minimizeVD)
import Numeric.LinearAlgebra qualified as LA
import Data.SRTree.Likelihoods
import Debug.Trace ( trace )

optimize :: Maybe Distribution -> Maybe Double -> Int -> Columns -> Column -> Maybe [Double] -> Fix SRTree -> (Fix SRTree, V.Vector Double, Int)
optimize mdist msErr niter xss ys mt0 tree = (tree', opt, its)
  where
    (tree', V.fromList -> t0') = floatConstsToParam tree
    t0          = case mt0 of
                    Nothing -> t0'
                    Just x  -> V.fromList x
    (opt, its)  = case mdist of
                    Nothing -> leastSquares niter xss ys tree' t0 
                    Just d  -> minimizeNLL d msErr niter xss ys tree' t0

leastSquares :: Int -> Columns -> Column -> Fix SRTree -> V.Vector Double -> (V.Vector Double, Int)
leastSquares niter xss ys tree t0
  | n == 0    = (t0, 0)
  | n > m     = (t0, 0)
  | otherwise = (toVec t_opt, iters path)
  where
    (t_opt, path) = nlFitting LevenbergMarquardtScaled 1e-6 1e-6 niter model jacob t0'
    toVec         = V.fromList . LA.toList
    iters         = fst . LA.size   
    t0'           = LA.fromList $ V.toList t0
    n             = LA.size t0'
    m             = LA.size ys
    model theta   = let theta' = V.fromList (LA.toList theta)
                     in subtract ys $ evalTree xss theta' LA.scalar tree
    jacob theta   = let theta' = V.fromList (LA.toList theta) 
                     in LA.fromColumns . snd $ gradParamsRev xss theta' LA.scalar tree

minimizeNLL :: Distribution -> Maybe Double -> Int -> Columns -> Column -> Fix SRTree -> V.Vector Double -> (V.Vector Double, Int)
minimizeNLL dist msErr niter xss ys tree t0
  | n == 0    = (t0, 0)
  | n > m     = (t0, 0)
  | otherwise = (toVec t_opt, iters path)
  where
    (t_opt, path) = minimizeVD VectorBFGS2 1e-6 niter 1e-3 1e-6 model jacob t0'

    toVec = V.fromList . LA.toList
    iters = fst . LA.size   
    t0'   = LA.fromList $ V.toList t0
    n     = LA.size t0'
    m     = LA.size ys
    model = nll dist msErr xss ys tree . V.fromList . LA.toList
    jacob = gradNLL dist msErr xss ys tree . V.fromList . LA.toList

minimizeGaussian :: Int -> Columns -> Column -> Fix SRTree -> V.Vector Double -> (V.Vector Double, Int)
minimizeGaussian = minimizeNLL Gaussian Nothing

minimizeBinomial :: Int -> Columns -> Column -> Fix SRTree -> V.Vector Double -> (V.Vector Double, Int)
minimizeBinomial = minimizeNLL Bernoulli Nothing
