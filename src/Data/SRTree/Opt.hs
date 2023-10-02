{-# language FlexibleInstances #-}
{-# language OverloadedStrings #-}
{-# language ImportQualifiedPost #-}
{-# language ViewPatterns #-}
module Data.SRTree.Opt
    ( optimize, sse, mse, rmse, Column, Columns, minimizeNLL, minimizeNLLWithFixedParam, minimizeGaussian, minimizeBinomial, minimizePoisson, nll, Distribution (..), gradNLL, fisherNLL )
    where

import Data.SRTree ( SRTree (..), Fix(..), floatConstsToParam )
import Data.SRTree.Eval ( evalTree )
import Data.SRTree.AD ( reverseModeUnique )
import Data.Vector.Storable qualified as VS
import Numeric.GSL.Fitting (FittingMethod (..), nlFitting)
import Numeric.GSL.Minimization (MinimizeMethodD (..), minimizeVD)
import Numeric.LinearAlgebra ( size, fromColumns )
import Data.SRTree.Likelihoods
import Debug.Trace ( trace )
import Data.SRTree.Print ( showExpr )

optimize :: Maybe Distribution -> Maybe Double -> Int -> Columns -> Column -> Maybe [Double] -> Fix SRTree -> (Fix SRTree, VS.Vector Double, Int)
optimize mDist mSErr nIter xss ys mTheta tree = (optTree, optTheta, steps)
  where
    (optTree, VS.fromList -> t0') = floatConstsToParam tree

    t0                = case mTheta of
                          Nothing -> t0'
                          Just x  -> VS.fromList x
    (optTheta, steps) = case mDist of
                          Nothing   -> leastSquares nIter xss ys optTree t0 
                          Just dist -> minimizeNLL dist mSErr nIter xss ys optTree t0

leastSquares :: Int -> Columns -> Column -> Fix SRTree -> VS.Vector Double -> (VS.Vector Double, Int)
leastSquares niter xss ys tree t0
  | niter == 0 = (t0, 0)
  | n == 0     = (t0, 0)
  | n > m      = (t0, 0)
  | otherwise  = (t_opt, iters path)
  where
    (t_opt, path) = nlFitting LevenbergMarquardtScaled 1e-20 1e-15 niter model jacob t0
    iters         = fst . size  
    n             = VS.length t0
    m             = VS.length ys
    model theta   = subtract ys $ evalTree xss theta VS.singleton tree
    jacob theta   = negate . fromColumns . snd $ reverseModeUnique xss theta VS.singleton tree

minimizeNLL :: Distribution -> Maybe Double -> Int -> Columns -> Column -> Fix SRTree -> VS.Vector Double -> (VS.Vector Double, Int)
minimizeNLL dist msErr niter xss ys tree t0
  | niter == 0 = (t0, 0)
  | n == 0     = (t0, 0)
  | n > m      = (t0, 0)
  | otherwise  = (t_opt, iters path)
  where
    (t_opt, path) = minimizeVD ConjugateFR 1e-20 niter 1e-1 1e-2 model jacob t0
    -- (t_opt, path) = minimizeVD ConjugateFR 1e-20 niter 1e-15 1e-3 model jacob t0

    iters = fst . size   
    n     = VS.length t0
    m     = VS.length ys
    model = nll dist msErr xss ys tree
    jacob = gradNLL dist msErr xss ys tree

minimizeNLLWithFixedParam :: Distribution -> Maybe Double -> Int -> Columns -> Column -> Fix SRTree -> Int -> VS.Vector Double -> (VS.Vector Double, Int)
minimizeNLLWithFixedParam dist msErr niter xss ys tree ix t0
  | niter == 0 = (t0, 0)
  | n == 0     = (t0, 0)
  | n > m      = (t0, 0)
  | otherwise  = (t_opt, iters path)
  where
    (t_opt, path) = minimizeVD ConjugateFR 1e-20 niter 1e-1 1e-2 model jacob t0
    --(t_opt, path) = minimizeVD ConjugateFR 1e-20 niter 1e-16 1e-3 model jacob t0

    iters   = fst . size   
    n       = VS.length t0
    m       = VS.length ys
    model   = nll dist msErr xss ys tree
    jacob t = gradNLL dist msErr xss ys tree t VS.// [(ix, 0.0)] 

minimizeGaussian :: Int -> Columns -> Column -> Fix SRTree -> VS.Vector Double -> (VS.Vector Double, Int)
minimizeGaussian = minimizeNLL Gaussian Nothing

minimizeBinomial :: Int -> Columns -> Column -> Fix SRTree -> VS.Vector Double -> (VS.Vector Double, Int)
minimizeBinomial = minimizeNLL Bernoulli Nothing

minimizePoisson :: Int -> Columns -> Column -> Fix SRTree -> VS.Vector Double -> (VS.Vector Double, Int)
minimizePoisson = minimizeNLL Poisson Nothing
