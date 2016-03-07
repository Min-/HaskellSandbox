{-#LANGUAGE OverloadedStrings#-}

{-
  Min Zhang
  5-29-2015

  import tables, mostly csv and txt, to list of lists

  all the text will be in Text type

   v.0.2.0 (July 24, 2015): add TableData type
   v.0.2.1 (Oct 9, 2015): add Type signature to all functions
-}

module MyTable
    ( importTable
    , smartTable
    , writeTable
    , column
    , cols
    , TableData
    , infinitp
    , infinitn
    , int
    , getInt
    , text
    , getText
    , double
    , getDouble
    , bool
    , getBool
--    , zipData
    , colType
    , select
    , getCol
    , getDoubleCol
    , getIntCol
    , getTextCol
  )
where

import qualified Data.Text.Lazy as T 
import qualified Data.Text.Lazy.IO as TextIO
import qualified Data.Char as C
import qualified Data.String as Str
import Control.Applicative
import qualified Data.List as L
import Control.Monad (fmap)
import Data.Ord (comparing)
import qualified Safe as S
import qualified Data.Foldable as F (all)

import qualified Data.ByteString.Lazy.Char8 as Bl
import qualified System.IO as IO

import MyText

--IO functions

-- if the table is windows format, then need to remove quotations too
importTable::(Str.IsString a, Eq a) => a -> FilePath -> IO [[T.Text]]
importTable delim input = 
  if delim == "\t" 
  then do map (tab . cleanWinNl . cleanQuote) . T.lines <$> TextIO.readFile input
  else if delim == ","
       then do  map (comma . cleanWinNl . cleanQuote) . T.lines <$> TextIO.readFile input
       else return []

-- use lazy bytestring to import first few lines (or first few bytes) of the file, not to overload the memory, just to text that type of table the import is. Therefore, to ease the use of importTable. Not sure if it's a good idea or not. Need to test it out.
smartTable :: FilePath -> IO [[T.Text]]
smartTable input = do
  inhandle <- IO.openFile input IO.ReadMode
  first1k <- Bl.hGet inhandle 1000 
  let t = Bl.count '\t' first1k
  let c = Bl.count ',' first1k
  let r = Bl.count '\r' first1k
  if t > c
  then do importTable "\t" input
  else do importTable "," input

writeTable :: FilePath -> [[T.Text]] -> IO ()
writeTable output table = do
  let result = T.unlines $ map untab table
  TextIO.writeFile output result


-- start from 1
column :: Int -> [a] -> [a]
column n = take 1 . drop (n - 1)

cols::Int -> [[a]] -> [a]
cols n = concat . map (column n)

-- TableData Type to put T.Text, Int, Double and Bool into the same type
data TableData
    = DataInt Int
    | DataText T.Text
    | DataDouble Double
    | DataBool Bool
    | Infp
    | Infn
       deriving (Show, Read, Eq, Ord)

infinitp::Double
infinitp = 1e100

infinitn::Double
infinitn = -1e100

--test1 :: [T.Text]
--test1 = ["1", "abc", "12.512"]             
--test2 = [DataInt 1, DataText "abc", DataDouble 12.512]

int :: T.Text -> TableData
int = DataInt . read . T.unpack

getInt :: TableData -> Int
getInt (DataInt a) = a

text :: T.Text -> TableData
text = DataText

getText :: TableData -> T.Text
getText (DataText a) = a

double :: T.Text -> TableData
double "Inf" = Infp
double "-Inf" = Infn
double x = (DataDouble . read . T.unpack) x

getDouble :: TableData -> Double
getDouble (DataDouble a) = a
getDouble Infp = infinitp
getDouble Infn = infinitn

bool :: T.Text -> TableData
bool = DataBool . read . T.unpack

getBool :: TableData -> Bool
getBool (DataBool a) = a

zipData :: [T.Text -> TableData] -> [T.Text] -> [TableData]
zipData f t = fmap (\(x,y) -> x y) $  zip f t

colType :: [T.Text -> TableData] -> [T.Text] -> [TableData]
colType = zipData

select :: (a -> b) -> (b -> Bool) -> [a] -> [a]
select col q = filter (q . col)

getCol :: (a -> b) -> [a] -> [b]
getCol col = map col 

getDoubleCol :: Int -> [TableData] -> Double
getDoubleCol n = getDouble . head . (column n)

getIntCol :: Int -> [TableData] -> Int
getIntCol n = getInt . head . (column n)

getTextCol :: Int -> [TableData] -> T.Text
getTextCol n = getText . head . (column n)


