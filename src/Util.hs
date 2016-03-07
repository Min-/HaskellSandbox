module Util
  ( groupCount
  , groupPct
  )
where

{-
Project name: Util library
Min Zhang
Date: Oct 9, 2015
Version: v.0.2.0
README: Add type signature for all the functions
-}


import qualified Data.List as L
import qualified Data.Text.Lazy as T
import Data.Ord (comparing)

-- utility functions
groupCount :: Ord a => [a] -> [(a, Int)]
groupCount x = zip (map head groupedOnes) (map length groupedOnes)
  where groupedOnes = groupDes x
       
groupPct :: (Ord a, Fractional b) => [a] -> [(a, b)] 
groupPct l = zip (map head groupedOnes) (map pct groupedOnes) 
  where groupedOnes = groupDes l
        totalLength = L.genericLength l
        pct x = L.genericLength x / totalLength

groupDes :: Ord a => [a] -> [[a]]
groupDes x = L.reverse . L.sortBy (comparing length) . L.group . L.reverse . L.sort $ x



{- not in use

binNumbers [] = []
binNumbers x = groupCount $ map toNearest x

toNearest x
  |0 > x = "Something's wrong"
  |0 <= x && x <= 1000 = "1k"
  |1000 < x && x <= 10000 = "10k"
  |10000 < x && x <= 50000 = "50k"
  |50000 < x && x <= 100000 = "100k"
  |100000 < x && x <= 500000 = "500k"
  |500000 < x && x <= 1000000 = "1M"
  |1000000 < x && x <= 10000000 = "10M"
  |10000000 < x = "over10M" 
-}

