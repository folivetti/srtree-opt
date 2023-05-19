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
import qualified Data.Vector.Storable as VS
import qualified Numeric.LinearAlgebra as LA
import Data.SRTree (SRTree (..), gradParams, gradParams, evalTree)
import Data.SRTree.Recursion (Fix(..))

type Columns = V.Vector Column
type Column  = LA.Vector Double

-- | Supported distributions for negative log-likelihood
data Distribution = Gaussian | Bernoulli 
    deriving (Show, Read, Enum, Bounded)

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
nll Bernoulli _ xss ys tree theta
  | notValid ys = error "For Bernoulli distribution the output must be either 0 or 1."
  | otherwise   = negate . LA.sumElements $ VS.zipWith (\a b -> if a == 0.0 then log (1 - b) else log b) ys yhat
  where
    yhat = evalTree xss theta LA.scalar tree
    notValid = VS.any (\x -> x /= 0 && x /= 1)
    

-- | Gradient of the negative log-likelihood
gradNLL :: Distribution -> Maybe Double -> Columns -> Column -> Fix SRTree -> V.Vector Double -> LA.Vector Double
gradNLL Gaussian msErr xss ys tree theta = LA.fromList [LA.sumElements (g * err) / sErr * sErr | g <-grad]
  where
    m            = fromIntegral $ LA.size ys
    n            = fromIntegral $ V.length theta
    (yhat, grad) = gradParams xss theta LA.scalar tree
    err          = yhat - ys
    sErr         = case msErr of
                     Nothing -> sqrt $ LA.sumElements (err ^ (2 :: Int)) / (m - n)
                     Just x  -> x
gradNLL Bernoulli _ xss ys tree theta
  | notValid ys = error "For Bernoulli distribution the output must be either 0 or 1."
  | otherwise   = LA.fromList [LA.sumElements $ (yhat - ys) * g / ((1 - yhat) * yhat) | g <- grad]
  where
    (yhat, grad) = gradParams xss theta LA.scalar tree
    notValid = VS.any (\x -> x /= 0 && x /= 1)
