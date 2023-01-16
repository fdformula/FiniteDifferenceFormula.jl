echo Testing 352 ...
julia Test_352.jl > a.tmp
diff a.tmp test_352.txt

echo Testing Wiki ...
julia Test_Wiki.jl > b.tmp
diff b.tmp test_wiki.txt

