module Data.SRTree.Likelihoods
  ( Distribution (..)
  , Column
  , Columns
  , sse
  , mse
  , rmse
  , nll
  , gradNLL
  , fisherNLL
  )
    where

import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS
import qualified Numeric.LinearAlgebra as LA
import Data.SRTree (SRTree (..), gradParamsRev, evalTree, deriveBy, floatConstsToParam)
import Data.SRTree.Recursion (Fix(..))
import Data.Maybe ( fromMaybe )

type Columns = V.Vector Column
type Column  = LA.Vector Double

-- | Supported distributions for negative log-likelihood
data Distribution = Gaussian | Bernoulli | Poisson
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

logit x = 1 / (1 + exp (-x))

-- | Negative log-likelihood
nll :: Distribution -> Maybe Double -> Columns -> Column -> Fix SRTree -> V.Vector Double -> Double
nll Gaussian msErr xss ys t theta = 0.5*ssr/(sErr*sErr) + 0.5*m*log (2*pi*sErr*sErr)
  where
    m   = fromIntegral $ LA.size ys
    ssr = sse xss ys t theta
    sErr = fromMaybe 1.0 msErr

-- | Bernoulli distribution of f(x; theta) is, given phi = 1 / (1 + exp (-f(x; theta))),
-- y log phi + (1-y) log (1 - phi), assuming y \in {0,1}
nll Bernoulli _ xss ys tree theta
  | notValid ys = error "For Bernoulli distribution the output must be either 0 or 1."
  | otherwise   = negate . LA.sumElements $ VS.zipWith (\a b -> if a == 0.0 then log (1 - b) else log b) ys yhat
  where
    yhat = logit $ evalTree xss theta LA.scalar tree
    notValid = VS.any (\x -> x /= 0 && x /= 1)

nll Poisson _ xss ys tree theta 
  | notValid ys = error "For Poisson distribution the output must be non-negative."
  | otherwise = negate . VS.sum $ ys * yhat - exp yhat
  where
    yhat = evalTree xss theta LA.scalar tree
    notValid = VS.any (<0)

-- | Gradient of the negative log-likelihood
gradNLL :: Distribution -> Maybe Double -> Columns -> Column -> Fix SRTree -> V.Vector Double -> LA.Vector Double
gradNLL Gaussian msErr xss ys tree theta = LA.fromList [LA.sumElements (g * err) / sErr * sErr | g <- grad]
  where
    (yhat, grad) = gradParamsRev xss theta LA.scalar tree
    err          = yhat - ys
    sErr         = fromMaybe 1.0 msErr

gradNLL Bernoulli _ xss ys tree theta
  | notValid ys = error "For Bernoulli distribution the output must be either 0 or 1."
  | otherwise   = LA.fromList [LA.sumElements $ (logit yhat - ys) * g | g <- grad]
  where
    (yhat, grad) = gradParamsRev xss theta LA.scalar tree
    notValid = VS.any (\x -> x /= 0 && x /= 1)

gradNLL Poisson _ xss ys tree theta
  | notValid ys = error "For Poisson distribution the output must be non-negative."
  | otherwise = LA.fromList [negate . LA.sumElements $ g * (ys - exp yhat) | g <- grad]
  where
    (yhat, grad) = gradParamsRev xss theta LA.scalar tree
    notValid = VS.any (<0)

-- | Fisher information of negative log-likelihood
fisherNLL :: Distribution -> Maybe Double -> Columns -> Column -> Fix SRTree -> V.Vector Double -> [Double]
fisherNLL d msErr xss ys tree theta = do 
  ix <- [0 .. p-1]
  let f'      = deriveBy True ix t'
      f''     = deriveBy True ix f'
      fvals'  = evalTree xss theta LA.scalar f'
      fvals'' = evalTree xss theta LA.scalar f''
  pure $ diagHess d sErr ys yhat fvals' fvals''
  where
    p       = V.length theta
    (t', _) = floatConstsToParam tree
    yhat    = evalTree xss theta LA.scalar t'
    ssr     = sse xss ys tree theta
    sErr    = case msErr of
                Nothing -> sqrt $ ssr / fromIntegral (LA.size ys - p)
                Just x  -> x

diagHess :: Distribution -> Double -> Column -> Column -> Column -> Column -> Double
diagHess Gaussian sErr ys yhat fvals' fvals'' = sum f_ii / (sErr ^ 2)
  where
    f_ii = LA.toList $ fvals'^2 - res*fvals''
    res  = ys - yhat

diagHess Bernoulli _ ys yhat fvals' fvals'' = sum . LA.toList $ (1 - phi)*phi*fvals'^2 + (phi - ys)*fvals''
  where phi = logit yhat

diagHess Poisson _ ys yhat fvals' fvals'' = sum . LA.toList $ phi*fvals'^2 + (phi - ys)*fvals''
  where phi = exp yhat

-- | Hessian of negative log-likelihood