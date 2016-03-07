{-#LANGUAGE OverloadedStrings#-}

{-
Project name: Lamina associated domain and find for domain boundaries and gene blocks
Min Zhang
Date: March 4, 2016
Version: v0.1.0
README: 
-}


-- Min Zhang
-- March 3, 2016

-- Find index of zero-crossing point

import qualified Safe as Safe

seq1 = [1,1,4,5,7,4,2,0,-1,-3,-9,-5,-2,4,12,15,17,-14,18,1,2,3,4]


import qualified Data.Text.Lazy as T
import qualified Data.Text.Lazy.IO as TextIO
import qualified Data.Char as C
import Control.Applicative
import qualified Data.List as L
import Control.Monad (fmap)
import Data.Ord (comparing)
import Data.Function (on)
import qualified Safe as Safe
import qualified Data.HashMap.Lazy as M
import qualified Data.Maybe as Maybe
import qualified Data.Foldable as F (all)
import Data.Traversable (sequenceA)
import qualified System.IO as IO
import System.Environment
import MyText
import MyTable
import Util

inputpath_sample = "data/GSE68260_Clone.5-5.1N.OE_LP150415.head.txt"

-- data structure: matrix with first row as colnames; first 3 columns are: chr start end
-- data contains NA
-- need to keep first 3 columns and do transformation functions (e.g. sum, mean etc) for the rest of columns


-- matrix transformation
-- matTrans n1 (number of columns to skip) n2 (number of rows to skip) d (delimiters) f (function for transformation) matrtix

matTrans n1 n2 d f mat = zipWith (\x y -> x ++ [y]) lables (trans mat3)
  where mat2 = drop n2 $ map (T.splitOn d) $ T.lines mat 
        mat3 = map (drop n1) mat2
        lables = map (take n1) mat2
        trans = map (readDouble . f . map (\x->if x == "NA" then 0 else toDouble x))

mean [] = 0
mean xs = sum xs / (L.genericLength xs)

zeroCrossing [] = []
zeroCrossing [x] = []
zeroCrossing  xs = 
  if product2 ys <= 0
  then 1 : zeroCrossing (drop 1 xs)  
  else 0 : zeroCrossing (drop 1 xs)
    where product2 x = (signum (fst $ head x)) * (signum (fst $ head $ drop 1 x))
          ys = zip xs [1..]
         
slideFunc f _ _ [] = [] 
slideFunc f windowSize slideSize xs
  | slideSize > windowSize = Safe.headNote "Slide size should be smaller than window size" []
  | length xs < windowSize = []
  | otherwise = (f $ take windowSize xs) : (slideFunc f windowSize slideSize $ drop slideSize xs)

slideSum = slideFunc sum

transformation = zeroCrossing . logBase 2 . mean

main = do
  input <- TextIO.readFile inputpath
  let res = matTrans 3 1 "\t" transformation input
  return res

