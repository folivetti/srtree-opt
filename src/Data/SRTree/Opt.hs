{-# language FlexibleInstances #-}
{-# language OverloadedStrings #-}
{-# language ImportQualifiedPost #-}
{-# language ViewPatterns #-}
module Data.SRTree.Opt
    ( optimize, sse, mse, rmse, Column, Columns, minimizeGaussian, minimizeBinomial, nll, Distribution (..), gradNLL, fisherNLL )
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
optimize mDist mSErr nIter xss ys mTheta tree = (optTree, optTheta, steps)
  where
    (optTree, V.fromList -> t0') = floatConstsToParam tree

    t0                = case mTheta of
                          Nothing -> t0'
                          Just x  -> V.fromList x
    (optTheta, steps) = case mDist of
                          Nothing   -> leastSquares nIter xss ys optTree t0 
                          Just dist -> minimizeNLL dist mSErr nIter xss ys optTree t0

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
    model = nll dist msErr xss ys tree . toVec
    jacob = gradNLL dist msErr xss ys tree . toVec

minimizeGaussian :: Int -> Columns -> Column -> Fix SRTree -> V.Vector Double -> (V.Vector Double, Int)
minimizeGaussian = minimizeNLL Gaussian Nothing

minimizeBinomial :: Int -> Columns -> Column -> Fix SRTree -> V.Vector Double -> (V.Vector Double, Int)
minimizeBinomial = minimizeNLL Bernoulli Nothing

minimizePoisson :: Int -> Columns -> Column -> Fix SRTree -> V.Vector Double -> (V.Vector Double, Int)
minimizePoisson = minimizeNLL Poisson Nothing