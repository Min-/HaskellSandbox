-- Min Zhang
-- March 2, 2016

-- LongestAlignment
-- BLAST algorithm

-- Find longest strech of matching
-- random generation sequences that match the probability
-- Not working yet



seq1 = [0, -1, -2, -1, -2, -1, -2, -1, 0, 1, 0, -1, -2, -3, -4] 

value = fst

pos = snd 

f [] (a, b, c, d) = (a, b, c, d)
f (x:xs) (min_, max_, minLoc, maxLoc)
  | value x < min_ && minLoc <= maxLoc = f xs (value x, max_, pos x, pos x)
  | value x < min_ = f xs (min_, max_, minLoc, maxLoc)
  | value x > max_ && (pos x) - minLoc + 1 > len = f xs (min_, value x, minLoc, pos x)
  | otherwise = f xs (min_, max_, minLoc, maxLoc)
  where len = maxLoc - minLoc + 1
