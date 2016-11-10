#!/bin/sh

javac Gotoh.java ScoringMatrix.java

for rep in $( seq 0 2 ); do
for init in 0 1; do
for int in 0 1; do
for mode in 1 2 3 4; do
	/usr/bin/time -a -o runtime.txt --format="$init $int $mode %U" java Gotoh --int $int --mode $mode --reinit $init >> /dev/null

done
done
done
done

Rscript plot.R

