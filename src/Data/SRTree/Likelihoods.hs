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
-- TODO: replace Data.Vector with Data.Vector.Storable
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

-- | Sum-of-square errors or Sum-of-square residues
sse :: Columns -> Column -> Fix SRTree -> V.Vector Double -> Double
sse xss ys tree theta = sum err
  where
    yhat  = evalTree xss theta LA.scalar tree
    err   = LA.toList $ (ys - yhat) ^ (2 :: Int)

-- | Mean squared errors
mse :: Columns -> Column -> Fix SRTree -> V.Vector Double -> Double
mse xss ys tree theta = sse xss ys tree theta / fromIntegral (LA.size ys)

-- | Root of the mean squared errors
rmse :: Columns -> Column -> Fix SRTree -> V.Vector Double -> Double
rmse xss ys tree = sqrt . mse xss ys tree

-- | logistic function
logistic x = 1 / (1 + exp (-x))
{-# inline logistic #-}

-- | get the standard error from a Maybe Double
-- if it is Nothing, estimate from the ssr, otherwise use the current value
-- For distributions other than Gaussian, it defaults to a constant 1
getSErr Gaussian est = fromMaybe est
getSErr _        _   = const 1
{-# inline getSErr #-}

-- negation of the sum of values in a vector
negSum = negate . LA.sumElements
{-# inline negSum #-}

-- | Negative log-likelihood
nll :: Distribution -> Maybe Double -> Columns -> Column -> Fix SRTree -> V.Vector Double -> Double

-- | Gaussian distribution
nll Gaussian msErr xss ys t theta = 0.5*(ssr/s2 + m*log (2*pi*s2))
  where
    m    = fromIntegral $ LA.size ys
    p    = fromIntegral $ V.length theta
    ssr  = sse xss ys t theta
    est  = sqrt $ ssr / (m - p)
    sErr = getSErr Gaussian est msErr
    s2   = sErr ^ 2

-- | Bernoulli distribution of f(x; theta) is, given phi = 1 / (1 + exp (-f(x; theta))),
-- y log phi + (1-y) log (1 - phi), assuming y \in {0,1}
nll Bernoulli _ xss ys tree theta
  | notValid ys = error "For Bernoulli distribution the output must be either 0 or 1."
  | otherwise   = negSum $ VS.zipWith (\a b -> if a == 0.0 then log (1 - b) else log b) ys yhat
  where
    yhat     = logistic $ evalTree xss theta LA.scalar tree
    notValid = VS.any (\x -> x /= 0 && x /= 1)

nll Poisson _ xss ys tree theta 
  | notValid ys = error "For Poisson distribution the output must be non-negative."
  | otherwise = negSum $ ys * yhat - exp yhat
  where
    yhat     = evalTree xss theta LA.scalar tree
    notValid = VS.any (<0)

-- | Gradient of the negative log-likelihood
gradNLL :: Distribution -> Maybe Double -> Columns -> Column -> Fix SRTree -> V.Vector Double -> LA.Vector Double
gradNLL Gaussian msErr xss ys tree theta = LA.fromList [LA.sumElements (g * err) / sErr * sErr | g <- grad]
  where
    m            = fromIntegral $ LA.size ys
    p            = fromIntegral $ V.length theta
    (yhat, grad) = gradParamsRev xss theta LA.scalar tree
    err          = yhat - ys
    ssr          = sse xss ys tree theta
    est          = sqrt $ ssr / (m - p)
    sErr         = getSErr Gaussian est msErr

gradNLL Bernoulli _ xss ys tree theta
  | notValid ys = error "For Bernoulli distribution the output must be either 0 or 1."
  | otherwise   = LA.fromList [LA.sumElements $ (logistic yhat - ys) * g | g <- grad]
  where
    (yhat, grad) = gradParamsRev xss theta LA.scalar tree
    notValid     = VS.any (\x -> x /= 0 && x /= 1)

gradNLL Poisson _ xss ys tree theta
  | notValid ys = error "For Poisson distribution the output must be non-negative."
  | otherwise   = LA.fromList [negate . LA.sumElements $ g * (ys - exp yhat) | g <- grad]
  where
    (yhat, grad) = gradParamsRev xss theta LA.scalar tree
    notValid     = VS.any (<0)

-- | Fisher information of negative log-likelihood
fisherNLL :: Distribution -> Maybe Double -> Columns -> Column -> Fix SRTree -> V.Vector Double -> [Double]
fisherNLL dist msErr xss ys tree theta = do 
  ix <- [0 .. p-1]
  let dtdix   = deriveBy True ix t'
      d2tdix2 = deriveBy True ix dtdix
      f'      = eval dtdix
      f''     = eval d2tdix2
  pure . (/sErr^2) . LA.sumElements $ phi' * f'^2 - res * f''
  where
    m       = LA.size ys
    p       = V.length theta
    (t', _) = floatConstsToParam tree
    eval    = evalTree xss theta LA.scalar

    yhat    = eval t'
    res     = ys - phi

    phi     = case dist of
                Gaussian  -> yhat
                Bernoulli -> logistic yhat
                Poisson   -> exp yhat
    phi'    = case dist of
                Gaussian  -> 1
                Bernoulli -> phi*(1 - phi)
                Poisson   -> phi

    ssr     = sse xss ys tree theta
    sErr    = getSErr dist est msErr
    est     = sqrt $ ssr / fromIntegral (m - p)
{-
diagHess :: Distribution -> Double -> Column -> Column -> Column -> Column -> Double
diagHess Gaussian sErr ys yhat fvals' fvals'' = sum f_ii / (sErr ^ 2)
  where
    f_ii = LA.toList $ fvals'^2 - res*fvals''
    res  = ys - yhat

diagHess Bernoulli _ ys yhat fvals' fvals'' = sum . LA.toList $ (1 - phi)*phi*fvals'^2 + (phi - ys)*fvals''
  where phi = logistic yhat

diagHess Poisson _ ys yhat fvals' fvals'' = sum . LA.toList $ phi*fvals'^2 + (phi - ys)*fvals''
  where phi = exp yhat
-}

-- | Hessian of negative log-likelihood