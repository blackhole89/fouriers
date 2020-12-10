import Data.List

pow 0 = [[]]
pow n = let sub = (pow (n-1)) in sub++(map (n-1:) sub)

unif 0 = [[]]
unif n = let sub = (unif (n-1)) in (map (-1:) sub)++(map (1:) sub)

chi cs xs = foldl (\p c -> (xs!!c)*p) 1 cs

dot n f g = (/(2^n)) $ foldl (+) 0 $ map (\x -> (f x)*(g x)) $ unif n

fourier n f = map (\cs -> ( cs, dot n f (chi cs) ) ) (pow n)
