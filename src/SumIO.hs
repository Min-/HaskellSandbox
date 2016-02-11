-- Feb 11, 2016
-- same IO program as cpp code, in Haskell
-- it seems like much less efficient than cpp code
-- compiled cpp exe is only 16k, compared to haskell exe 1.8MB
-- but the good thing is: all the data is preserved, you can do more computation if necessary, but cpp code is mutated and raw data is lost.


-- both codes are easily broken if the input age is not a number

-- So in long term, Haskell code is definitely safer and might be more flexible and efficient.

import Control.Applicative
import qualified Safe 

main = do
  allages <- getAge
  putStrLn (show $ sum allages)
  putStrLn (show $ getAverage allages)

getAge = do
  putStrLn "put a number, -1 to terminate. "
  age <- checkNumber <$> getLine
  if age == -1
  then return []
  else do 
          list <- getAge
          return (age : list)

checkNumber :: String -> Int
checkNumber n = Safe.readNote "Not a number" n
-- or Safe.readDef 0 n

getAverage xs = sum xs `div` length xs

-- IO interactive program template from Erik's online class

--interactive prelude evaluate terminate success next word = do
--  prelude
--  xs <- getLine
--  if evaluate xs word
--  then terminate xs
--  else do success word xs
--          interactive prelude evaluate terminate success next (next word)
--  

  
  
