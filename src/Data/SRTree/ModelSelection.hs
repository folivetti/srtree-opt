module Data.SRTree.ModelSelection where

import qualified Data.Vector.Storable as VS
import Data.SRTree
import Data.SRTree.Recursion ( cata ) 
import Data.SRTree.Opt
import qualified Numeric.LinearAlgebra as LA

-- | Bayesian information criterion
bic :: Distribution -> Maybe Double -> Columns -> Column -> VS.Vector Double -> Fix SRTree -> Double
bic dist mSErr xss ys theta tree = (p + 1) * log n + 2 * nll dist mSErr xss ys tree theta
  where
    p = fromIntegral $ VS.length theta
    n = fromIntegral $ VS.length ys
{-# INLINE bic #-}

-- | Akaike information criterion
aic :: Distribution -> Maybe Double -> Columns -> Column -> VS.Vector Double -> Fix SRTree -> Double
aic dist mSErr xss ys theta tree = 2 * (p + 1) + 2 * nll dist mSErr xss ys tree theta
  where
    p = fromIntegral $ VS.length theta
    n = fromIntegral $ VS.length ys
{-# INLINE aic #-}

evidence :: Distribution -> Maybe Double -> Columns -> Column -> VS.Vector Double -> Fix SRTree -> Double
evidence dist mSErr xss ys theta tree = (1 - b) * nll dist mSErr xss ys tree theta - p / 2 * log b
  where
    p = fromIntegral $ VS.length theta
    n = fromIntegral $ VS.length ys
    b = 1 / sqrt n
{-# INLINE evidence #-}

-- | MDL as described in 
-- Bartlett, Deaglan J., Harry Desmond, and Pedro G. Ferreira. "Exhaustive symbolic regression." IEEE Transactions on Evolutionary Computation (2023).
mdl :: Distribution -> Maybe Double -> Columns -> Column -> VS.Vector Double -> Fix SRTree -> Double
mdl dist mSErr xss ys theta tree = nll' dist mSErr xss ys theta' tree
                                 + logFunctional tree
                                 + logParameters dist mSErr xss ys theta tree
  where
    fisher = fisherNLL dist mSErr xss ys tree theta
    theta' = VS.fromList $ map (\(t,f) -> if isSignificant t f then t else 0.0) $ zip (VS.toList theta) fisher
    isSignificant v f = abs (v / sqrt(12 / f) ) >= 1
{-# INLINE mdl #-}

mdlLatt :: Distribution -> Maybe Double -> Columns -> Column -> VS.Vector Double -> Fix SRTree -> Double
mdlLatt dist mSErr xss ys theta tree = nll' dist mSErr xss ys theta' tree
                                     + logFunctional tree
                                     + logParametersLatt dist mSErr xss ys theta tree
  where
    fisher = fisherNLL dist mSErr xss ys tree theta
    theta' = VS.fromList $ map (\(t,f) -> if isSignificant t f then t else 0.0) $ zip (VS.toList theta) fisher
    isSignificant v f = abs (v / sqrt(12 / f) ) >= 1
{-# INLINE mdlLatt #-}


-- | same as `mdl` but weighting the functional structure by frequency calculated using a wiki information of
-- physics and engineering functions
mdlFreq :: Distribution -> Maybe Double -> Columns -> Column -> VS.Vector Double -> Fix SRTree -> Double
mdlFreq dist mSErr xss ys theta tree = nll' dist mSErr xss ys theta tree
                                     + logFunctionalFreq tree
                                     + logParameters dist mSErr xss ys theta tree
{-# INLINE mdlFreq #-}

logFunctional :: Fix SRTree -> Double
logFunctional tree = countNodes tree * log (countUniqueTokens tree') 
                   + foldr (\c acc -> log (abs c) + acc) 0 consts 
                   + log(2) * numberOfConsts
  where
    tree'          = fst $ floatConstsToParam tree
    consts         = getIntConsts tree
    numberOfConsts = fromIntegral $ length consts
    signs          = sum [1 | a <- getIntConsts tree, a < 0] -- TODO: will we use that?
{-# INLINE logFunctional #-}

logFunctionalFreq  :: Fix SRTree -> Double
logFunctionalFreq tree = treeToNat tree' 
                       + foldr (\c acc -> log (abs c) + acc) 0 consts  
                       + countVarNodes tree * log (numberOfVars tree)
  where
    tree'  = fst $ floatConstsToParam tree
    consts = getIntConsts tree
{-# INLINE logFunctionalFreq #-}

logParameters :: Distribution -> Maybe Double -> Columns -> Column -> VS.Vector Double -> Fix SRTree -> Double
logParameters dist mSErr xss ys theta tree = -(p / 2) * log 3 + 0.5 * logFisher + logTheta
  where
    -- p      = fromIntegral $ VS.length theta
    fisher = fisherNLL dist mSErr xss ys tree theta

    (logTheta, logFisher, p) = foldr addIfSignificant (0, 0, 0)
                             $ zip (VS.toList theta) fisher

    addIfSignificant (v, f) (acc_v, acc_f, acc_p)
       | isSignificant v f = (acc_v + log (abs v), acc_f + log f, acc_p + 1)
       | otherwise         = (acc_v, acc_f, acc_p)

    isSignificant v f = abs (v / sqrt(12 / f) ) >= 1

logParametersLatt :: Distribution -> Maybe Double -> Columns -> Column -> VS.Vector Double -> Fix SRTree -> Double
logParametersLatt dist mSErr xss ys theta tree = 0.5 * p * (1 - log 3) + 0.5 * log detFisher
  where
    fisher = fisherNLL dist mSErr xss ys tree theta
    detFisher = LA.det $ LA.fromLists $ hessianNLL dist mSErr xss ys tree theta

    (logTheta, logFisher, p) = foldr addIfSignificant (0, 0, 0)
                             $ zip (VS.toList theta) fisher

    addIfSignificant (v, f) (acc_v, acc_f, acc_p)
       | isSignificant v f = (acc_v + log (abs v), acc_f + log f, acc_p + 1)
       | otherwise         = (acc_v, acc_f, acc_p)

    isSignificant v f = abs (v / sqrt(12 / f) ) >= 1

-- flipped version of nll
nll' :: Distribution -> Maybe Double -> Columns -> Column -> VS.Vector Double -> Fix SRTree -> Double
nll' dist mSErr xss ys theta tree = nll dist mSErr xss ys tree theta
{-# INLINE nll' #-}

treeToNat :: Fix SRTree -> Double
treeToNat = cata alg
  where
    alg (Uni f t)    = funToNat f + t
    alg (Bin op l r) = opToNat op + l + r
    alg _            = 0.6610799229372109

    opToNat :: Op -> Double
    opToNat Add = 2.500842464597881
    opToNat Sub = 2.500842464597881
    opToNat Mul = 1.720356134912558
    opToNat Div = 2.60436883851265
    opToNat Power = 2.527957363394847

    funToNat :: Function -> Double
    funToNat Sqrt = 4.780867285331753
    funToNat Log  = 4.765599813200964
    funToNat Exp  = 4.788589331425663
    funToNat Abs  = 6.352564869783006
    funToNat Sin  = 5.9848400896576885
    funToNat Cos  = 5.474014465891698
    funToNat Sinh = 8.038963823353235
    funToNat Cosh = 8.262107374667444
    funToNat Tanh = 7.85664226655928
    funToNat Tan  = 8.262107374667444
    funToNat _    = 8.262107374667444
    --funToNat Factorial = 7.702491586732021
{-# INLINE treeToNat #-}
