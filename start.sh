#!/bin/bash
for i in {5..95..5}; do
	mkdir ../../data/test/$i
	if [ ! -f ../../data/test/$i/reynoldsModelList.txt ]; then
		./testbed ../../../../data/TIMIT/timit-reynolds-train.crp - - ../../data/metagmm_dd.ini 25.0 17 $i -debug.debugDir=../../data/test/$i
	fi
	./testbed ../../../../data/TIMIT/timit-reynolds-test.crp - - ../../data/metagmm_dd.ini 25.0 18 ../../data/test/$i/reynoldsModelList.txt $i -debug.debugDir=../../data/test/$i
done

for i in {5..95..5}; do
        mkdir ../../data/test/$i/std
	if [ ! -f ../../data/test/$i/std/reynoldsModelList.txt ]; then
	        ./testbed ../../../../data/TIMIT/timit-reynolds-train.crp - - ../../data/metagmm_std.ini 25.0 17 $i -debug.debugDir=../../data/test/$i/std
	fi
        ./testbed ../../../../data/TIMIT/timit-reynolds-test.crp - - ../../data/metagmm_std.ini 25.0 18 ../../data/test/$i/std/reynoldsModelList.txt $i -debug.debugDir=../../data/test/$i/std
done

for i in {5..95..5}; do
        mkdir ../../data/test/$i/fix
       if [ ! -f ../../data/test/$i/fix/reynoldsModelList.txt ]; then
               ./testbed ../../../../data/TIMIT/timit-reynolds-train.crp - - ../../data/metagmm_fix.ini 25.0 17 $i -debug.debugDir=../../data/test/$i/fix
       fi
        ./testbed ../../../../data/TIMIT/timit-reynolds-test.crp - - ../../data/metagmm_fix.ini 25.0 18 ../../data/test/$i/fix/reynoldsModelList.txt $i -debug.debugDir=../../data/test/$i/fix
done


#if [ ! -f ../../data/test/80/reynoldsModelList.txt ]; then
#       ./testbed ../../../../data/TIMIT/timit-reynolds-train.crp - - ../../data/metagmm_dd.ini 25.0 17 80 -debug.debugDir=../../data/test/80
#fi

#./testbed ../../../../data/TIMIT/timit-reynolds-test.crp - - ../../data/metagmm_dd.ini 25.0 18 ../../data/test/80/reynoldsModelList.txt 80 -debug.debugDir=../../data/test/80
