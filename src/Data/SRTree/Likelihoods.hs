module Data.SRTree.Likelihoods 
  ( Distribution (..)
  , Column
  , Columns
  , sse
  , mse
  , rmse
  , nll
  , gradNLL) 
    where

import qualified Data.Vector as V
import qualified Numeric.LinearAlgebra as LA
import Data.SRTree (SRTree (..), gradParams, floatConstsToParam, gradParams, evalTree)
import Data.SRTree.Recursion (Fix(..), cata)

type Columns = V.Vector Column
type Column  = LA.Vector Double

-- | Supported distributions for negative log-likelihood
data Distribution = Gaussian | Binomial deriving (Show, Read, Enum, Bounded)

sse :: Columns -> Column -> Fix SRTree -> V.Vector Double -> Double
sse xss ys tree theta = sum err
  where
    yhat  = evalTree xss theta LA.scalar tree
    err   = LA.toList $ (ys - yhat) ^ (2 :: Int)

mse :: Columns -> Column -> Fix SRTree -> V.Vector Double -> Double
mse xss ys tree theta = sse xss ys tree theta / fromIntegral (LA.size ys)

rmse :: Columns -> Column -> Fix SRTree -> V.Vector Double -> Double
rmse xss ys tree = sqrt . mse xss ys tree

-- | Negative log-likelihood
nll :: Distribution -> Maybe Double -> Columns -> Column -> Fix SRTree -> V.Vector Double -> Double
nll Gaussian msErr xss ys t theta = 0.5*ssr/(sErr*sErr) + 0.5*m*log(2*pi*sErr*sErr)
  where
    m   = fromIntegral $ LA.size ys
    n   = fromIntegral $ V.length theta
    ssr = sse xss ys t theta
    sErr = case msErr of
           Nothing -> sqrt $ ssr / (m - n)
           Just x  -> x

-- | Gradient of the negative log-likelihood
gradNLL :: Distribution -> Maybe Double -> Columns -> Column -> Fix SRTree -> V.Vector Double -> LA.Vector Double
gradNLL Gaussian msErr xss ys tree theta = LA.fromList [LA.sumElements (g * err) / sErr^2 | g <-grad]
  where
    m            = fromIntegral $ LA.size ys
    n            = fromIntegral $ V.length theta
    (yhat, grad) = gradParams xss theta LA.scalar tree
    err          = yhat - ys
    sErr         = case msErr of
                     Nothing -> sqrt $ LA.sumElements (err ^ 2) / (m - n)
                     Just x  -> x
