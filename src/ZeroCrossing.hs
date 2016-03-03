-- Min Zhang
-- March 3, 2016

-- Find index of zero-crossing point

seq1 = [1,1,4,5,7,4,2,0,-1,-3,-9,-5,-2,4,12,15,17,-14,18,1,2,3,4]

zeroCrossing [] = []
zeroCrossing [x] = []
zeroCrossing  xs = 
  if product2 ys <= 0
  then 1 : zeroCrossing (drop 1 xs)  
  else 0 : zeroCrossing (drop 1 xs)
    where product2 x = (fst $ head x) * (fst $ head $ drop 1 x)
          ys = zip xs [1..]
          index = snd . head
          
      
