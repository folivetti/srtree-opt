module Data.SRTree.Likelihoods
  ( Distribution (..)
  , Column
  , Columns
  , sse
  , mse
  , rmse
  , nll
  , predict
  , gradNLL
  , fisherNLL
  , getSErr
  , hessianNLL
  )
    where
-- TODO: replace Data.Vector with Data.Vector.Storable
-- import qualified Data.Vector as V
import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS
import qualified Numeric.LinearAlgebra as LA
import Data.SRTree (SRTree (..), Fix(..), floatConstsToParam)
import Data.SRTree.AD -- ( reverseModeUnique )
import Data.SRTree.Eval ( evalTree )
import Data.SRTree.Derivative ( deriveByParam )
import Data.Maybe ( fromMaybe )

import Debug.Trace ( trace )

type Columns = V.Vector Column
type Column  = LA.Vector Double

-- | Supported distributions for negative log-likelihood
data Distribution = Gaussian | Bernoulli | Poisson
    deriving (Show, Read, Enum, Bounded)

-- | Sum-of-square errors or Sum-of-square residues
sse :: Columns -> Column -> Fix SRTree -> VS.Vector Double -> Double
sse xss ys tree theta = err
  where
    yhat  = evalTree xss theta VS.singleton tree
    err   = VS.sum $ (ys - yhat) ^ (2 :: Int)

-- | Mean squared errors
mse :: Columns -> Column -> Fix SRTree -> VS.Vector Double -> Double
mse xss ys tree theta = sse xss ys tree theta / fromIntegral (VS.length ys)

-- | Root of the mean squared errors
rmse :: Columns -> Column -> Fix SRTree -> VS.Vector Double -> Double
rmse xss ys tree = sqrt . mse xss ys tree

-- | logistic function
logistic :: Floating a => a -> a
logistic x = 1 / (1 + exp (-x))
{-# inline logistic #-}

-- | get the standard error from a Maybe Double
-- if it is Nothing, estimate from the ssr, otherwise use the current value
-- For distributions other than Gaussian, it defaults to a constant 1
getSErr :: Num a => Distribution -> a -> Maybe a -> a
getSErr Gaussian est = fromMaybe est
getSErr _        _   = const 1
{-# inline getSErr #-}

-- negation of the sum of values in a vector
negSum :: VS.Vector Double -> Double
negSum = negate . VS.sum
{-# inline negSum #-}

-- | Negative log-likelihood
nll :: Distribution -> Maybe Double -> Columns -> Column -> Fix SRTree -> VS.Vector Double -> Double

-- | Gaussian distribution
nll Gaussian msErr xss ys t theta = 0.5*(ssr/s2 + m*log (2*pi*s2))
  where
    m    = fromIntegral $ VS.length ys
    p    = fromIntegral $ VS.length theta
    ssr  = sse xss ys t theta
    est  = sqrt $ ssr / (m - p)
    sErr = getSErr Gaussian est msErr
    s2   = sErr ^ 2

-- | Bernoulli distribution of f(x; theta) is, given phi = 1 / (1 + exp (-f(x; theta))),
-- y log phi + (1-y) log (1 - phi), assuming y \in {0,1}
nll Bernoulli _ xss ys tree theta
  | notValid ys = error "For Bernoulli distribution the output must be either 0 or 1."
  | otherwise   = negSum $ VS.zipWith (\a b -> a*b - log (1 + exp b)) ys yhat
  -- | otherwise   = negSum $ VS.zipWith (\a b -> if a == 0.0 then log (1 - b) else log b) ys yhat
  where
    yhat     = evalTree xss theta VS.singleton tree
    -- yhat     = logistic $ evalTree xss theta VS.singleton tree
    notValid = VS.any (\x -> x /= 0 && x /= 1)

nll Poisson _ xss ys tree theta 
  | notValid ys = error "For Poisson distribution the output must be non-negative."
  | VS.any isNaN yhat = error $ "NaN predictions " <> show theta
  | otherwise   = negSum $ ys * yhat - ys * log ys - exp yhat
  where
    yhat     = evalTree xss theta VS.singleton tree
    notValid = VS.any (<0)

-- | Prediction for different distributions
predict :: Distribution -> Fix SRTree -> VS.Vector Double -> Columns -> LA.Vector Double
predict Gaussian  tree theta xss = evalTree xss theta VS.singleton tree
predict Bernoulli tree theta xss = logistic $ evalTree xss theta VS.singleton tree
predict Poisson   tree theta xss = exp $ evalTree xss theta VS.singleton tree

-- | Gradient of the negative log-likelihood
gradNLL :: Distribution -> Maybe Double -> Columns -> Column -> Fix SRTree -> VS.Vector Double -> VS.Vector Double
gradNLL Gaussian msErr xss ys tree theta = VS.fromList [VS.sum (g * err) / sErr * sErr | g <- grad]
  where
    m            = fromIntegral $ VS.length ys
    p            = fromIntegral $ VS.length theta
    -- (yhat, grad) = reverseModeUnique xss theta VS.singleton tree
    yhat         = predict Gaussian tree theta xss
    grad         = forwardMode xss theta VS.singleton tree
    err          = yhat - ys
    ssr          = sse xss ys tree theta
    est          = sqrt $ ssr / (m - p)
    sErr         = getSErr Gaussian est msErr

gradNLL Bernoulli _ xss ys tree theta
  | notValid ys = error "For Bernoulli distribution the output must be either 0 or 1."
  | otherwise   = VS.fromList [VS.sum $ (yhat - ys) * g | g <- grad']
  where
    m            = fromIntegral $ LA.size yhat
    yhat         = predict Bernoulli tree theta xss
    grad         = forwardMode xss theta VS.singleton tree
    notValid     = VS.any (\x -> x /= 0 && x /= 1)
    grad'        = map (LA.cmap nanTo0) grad
    nanTo0 x     = if isNaN x then 0 else x

gradNLL Poisson _ xss ys tree theta
  | notValid ys = error "For Poisson distribution the output must be non-negative."
  | any (VS.any isNaN) grad = error $ "NaN gradient " <> show grad
  | otherwise   = VS.fromList [VS.sum $ g * (yhat - ys) | g <- grad]
  where
    yhat         = predict Poisson tree theta xss
    grad         = forwardMode xss theta VS.singleton tree
    notValid     = VS.any (<0)

-- | Fisher information of negative log-likelihood
fisherNLL :: Distribution -> Maybe Double -> Columns -> Column -> Fix SRTree -> VS.Vector Double -> [Double]
fisherNLL dist msErr xss ys tree theta = do 
  ix <- [0 .. p-1]
  let dtdix   = deriveByParam ix t'
      d2tdix2 = deriveByParam ix dtdix
      f'      = eval dtdix
      f''     = eval d2tdix2
  pure . (/sErr^2) . VS.sum $ phi' * f'^2 - res * f''
  where
    m    = VS.length ys
    p    = VS.length theta
    t'   = fst $ floatConstsToParam tree
    eval = evalTree xss theta VS.singleton
    ssr  = sse xss ys tree theta
    sErr = getSErr dist est msErr
    est  = sqrt $ ssr / fromIntegral (m - p)
    yhat = eval t'
    res  = ys - phi

    (phi, phi') = case dist of
                    Gaussian  -> (yhat, 1)
                    Bernoulli -> (logistic yhat, phi*(1 - phi))
                    Poisson   -> (exp yhat, phi)

-- | Hessian of negative log-likelihood
--
-- Note, though the Fisher is just the diagonal of the return of this function
-- it is better to keep them as different functions for efficiency
hessianNLL :: Distribution -> Maybe Double -> Columns -> Column -> Fix SRTree -> VS.Vector Double -> [[Double]]
hessianNLL dist msErr xss ys tree theta = do 
  ix <- [0 .. p-1]
  let dtdix = deriveByParam ix t'
      fx    = eval dtdix
      rs    = do iy <- [0 .. p-1]
                 let d2tdixy = deriveByParam iy dtdix
                     dtdiy   = deriveByParam iy t'
                     fy      = eval dtdiy
                     fxy     = eval d2tdixy
                 pure . (/sErr^2) . VS.sum $ phi' * fx * fy - res * fxy
  pure rs
  where
    m    = VS.length ys
    p    = VS.length theta
    t'   = fst $ floatConstsToParam tree
    eval = evalTree xss theta VS.singleton
    ssr  = sse xss ys tree theta
    sErr = getSErr dist est msErr
    est  = sqrt $ ssr / fromIntegral (m - p)
    yhat = eval t'
    res  = ys - phi

    (phi, phi') = case dist of
                    Gaussian  -> (yhat, 1)
                    Bernoulli -> (logistic yhat, phi*(1 - phi))
                    Poisson   -> (exp yhat, phi)
