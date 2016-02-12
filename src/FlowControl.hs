-- for loops in haskell

import Control.Applicative

-- main1 print 1 to 9
-- in general, we can replace for loop with map
main1 = do
  mapM_ (\x -> print x) [1..9]

nTimes :: Int -> IO () -> IO ()
nTimes 0 do_this = return ()
nTimes n do_this = do
    do_this
    nTimes (n-1) do_this

main2 = do
  mapM_ prettyPrint [1..30]
  where prettyPrint x = putStrLn $ (show x) ++ " ---- " ++ (show $ investmentReturn 10000 0.01 x) 

investmentReturn principle rate day = principle * (1+rate)^day

main3 = do
  x <- (\n -> read n::Int) <$> getLine
  dowhile x

-- similar to nTimes
dowhile x = do
  print x
  if x >= 40
  then return ()
  else dowhile (x + 1)

getNumber = getLine >>= return . (\x->read x::Int)

-- switch or case
main = do
  x <- getNumber
  case x of
    16 -> putStrLn "buy beer"
    21 -> putStrLn "lottery!!"
    _ -> print "default"
