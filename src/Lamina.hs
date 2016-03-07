{-#LANGUAGE OverloadedStrings#-}
{-#LANGUAGE FlexibleContexts#-}

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

-- ## Samples
seq1 = [1,1,4,5,7,4,2,0,-1,-3,-9,-5,-2,4,12,15,17,-14,18,1,2,3,4]

-- sequencing data
inputpath_sample = "/Users/minzhang/Documents/data/P55_hiC_looping/data/GSE68260_Clone.5-5.1N.OE_LP150415.head.txt"
inputpath = "/Users/minzhang/Documents/data/P55_hiC_looping/data/GSE68260_Clone.5-5.1N.OE_LP150415.txt"

-- array data
-- TODO: note that the ref is using hg18, liftOver the final wig, and do wigToBigWig
--input_array_ref = "/Users/minzhang/Documents/data/P55_hiC_looping/data/GPL10559-18779.txt"
--output_intermediate_ref = "/Users/minzhang/Documents/data/P55_hiC_looping/data/GPL10559-18779.clean.txt"

--input_array_data = "/Users/minzhang/Documents/data/P55_hiC_looping/data/GSM1700333-6918.txt"
--output_annotated_array = "/Users/minzhang/Documents/data/P55_hiC_looping/data/GSM1700333-6918.array.txt"

-- ## basic functions
mean [] = 0
mean xs = sum xs / (L.genericLength xs)

zeroCrossing [] = []
zeroCrossing [x] = []
zeroCrossing  xs 
  | product2 ys && fstSig ys >= 0 = "+" : zeroCrossing (drop 1 xs)  
  | product2 ys && fstSig ys < 0  = "-" : zeroCrossing (drop 1 xs)
  | otherwise                     = "*" : zeroCrossing (drop 1 xs)
    where product2 x = fstSig x == sndSig x
          ys = zip xs [1..]
          fstSig x = signum (fst $ head x)
          sndSig x = signum (fst $ head $ drop 1 x)
         
slideFunc f _ _ [] = [] 
slideFunc f windowSize slideSize xs
  | slideSize > windowSize = Safe.headNote "Slide size should be smaller than window size" []
  | length xs < windowSize = []
  | otherwise = (f $ take windowSize xs) : (slideFunc f windowSize slideSize $ drop slideSize xs)

slideSum = slideFunc sum

transformation = mean

-- data structure: matrix with first row as colnames; first 3 columns are: chr start end
-- data contains NA, NAs are turned into zeros
-- need to keep first 3 columns and do transformation functions (e.g. sum, mean etc) for the rest of columns

-- matrix transformation
-- matTrans n1 (number of columns to skip) n2 (number of rows to skip) d (delimiters) f (function for transformation) matrtix
matTrans :: Int -> Int -> T.Text -> ([Double] -> b) -> T.Text -> [([T.Text], b)]
matTrans n1 n2 d f mat = zip labels (trans mat3)
  where mat2 = drop n2 $ map (T.splitOn d) $ T.lines mat 
        mat3 = map (drop n1) mat2
        labels = map (take n1) mat2
        trans = map (f . map (\x->if x == "NA" then 0 else toDouble x))

-- specific example: get a chromosome data and plot in UCSC
-- there are special chr names that invalid in this data, needs to be cleaned up
pair2Wig :: Foldable t => [([T.Text], Double)] -> t T.Text -> T.Text
pair2Wig xs filterChr = combineEachChr . filterList filterChr . groupList . sortList $ xs
  where
    combineEachChr xs = T.unlines $ zipWith (\x y -> T.init $ T.unlines [chrHeader x, y]) (cleanEachChr xs) (cleanEachData xs)
    cleanEachChr xs = map (head . fst . head) xs 
    -- take grouped positions (grouped by chr), each element is a tuple, ([chr, start, end], value)
    cleanEachData xs = map (T.init . T.unlines . map (\x-> T.intercalate "\t" [head $ drop 1 $ fst x, (readDouble $ snd x)])) xs
    --                                                                        [start site (drop end site), value] -- for array, some probes have potential same start sites
    filterList chrs = filter (\x-> ( not . (\y -> y `elem` chrs) . head .fst . head) x) -- remove invalid chromosomes
    groupList = L.groupBy ((==) `on` (head . fst)) -- group by Chr
    sortList = L.sortBy (comparing (head . fst))  -- sort by Chr

wigHeader = "track type=wiggle_0 name=\"variableStep\" description=\"variableStep format\" visibility=full autoScale=off viewLimits=0.0:25.0 color=50,150,255 priority=10\n"
chrHeader chr = T.concat ["variableStep chrom=", chr]

-- ## main function to produce results
-- LamID-seq single cell matrix data to wig.
matrix2wig = do
  input <- TextIO.readFile inputpath
  let outputpath = "/Users/minzhang/Documents/data/P55_hiC_looping/data/GSE68260_Clone.5-5.1N.OE_LP150415.wig" 
  let res = matTrans 3 1 "\t" transformation input
  TextIO.writeFile outputpath wigHeader
  TextIO.appendFile outputpath (pair2Wig res ["chr9:22", "chr22:9"])

-- March 7, 2016
-- Turn array data into wig
-- generate a intermediate file that only keep probeID, chr, start (since most of the probes are short enough)
arrayRef2CleanFile input_array_ref output_intermediate_ref = do
  input <- map (T.splitOn "\t") . dropWhile (\x->T.head x == '#' || T.take 2 x == "ID") . T.lines . T.replace "\r" "" <$> TextIO.readFile input_array_ref  
  let res = map (\x->[x!!0, x!!1, x!!3, x!!4]) input
  writeTable output_intermediate_ref res

-- need to update arrayRef2CleanFile first, Only need to do it once
arrayRef2Map output_intermediate_ref = do
  input <- smartTable output_intermediate_ref
  return $ M.fromList $ map (\x-> (x!!0, [x!!1, x!!2, x!!3])) input

-- annotate array data
annotateArray input_ref input_array_data output_annotated_array = do
  array <- map (T.splitOn "\t") . dropWhile (\x->T.head x == '#' || T.take 2 x == "ID") . T.lines . T.replace "\r" "" <$> TextIO.readFile input_array_data
  ref <- arrayRef2Map input_ref
  let res = filter (\x->head x /= "NA") $ map (\x-> (M.lookupDefault ["NA"] (head x) ref) ++ (tail x)) array
  writeTable output_annotated_array res

array2Wig input_sorted_bed output_wig = do
 -- let inputpath = "/Users/minzhang/Documents/data/P55_hiC_looping/data/GSM1700333-6918.array.hg19.sort.bed" -- liftOver from hg18 to hg19
 -- let outputpath = "/Users/minzhang/Documents/data/P55_hiC_looping/data/GSM1700333-6918.hg19.wig"
  valuePair <- map (\x->(take 3 x, toDouble $ last x)) <$> smartTable input_sorted_bed
  TextIO.writeFile output_wig wigHeader
  TextIO.appendFile output_wig (pair2Wig valuePair []) -- good, re-use code

-- compile array2Wig
main = do
  TextIO.putStrLn "Lamina input_ref input_array_data output_annotated_array" 
  [mode, refpath, input, output] <- take 4 <$> getArgs
  case mode of 
    "array2bed" -> do annotateArray refpath input output
    "bed2wig" -> do array2Wig input output
