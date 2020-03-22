#!/bin/bash
g++ -Wall -pedantic -fopenmp -Ofast -std=c++11 -o mhr mhr.cpp
./mhr 14.612 43000 0.01 < testData/mhr_20_10_5.txt
./mhr 22.219 17000000 0.4 < testData/mhr_30_10_10.txt
./mhr 28.245 76000000 1.8 < testData/mhr_30_10_15.txt	
./mhr 24.704 210000000 5.5 < testData/mhr_34_10_15.txt
./mhr 29.326 291000000 7.3 < testData/mhr_34_10_17.txt
./mhr 55.436 6000000000 178 < testData/mhr_37_15_17.txt
./mhr 56.816 14000000000 428 < testData/mhr_37_15_18.txt
./mhr 56.738 14000000000 410 < testData/mhr_38_15_18.txt	
./mhr 74.258 9000000000 276 < testData/mhr_38_20_15.txt	
./mhr 87.362 43000000000 1342 < testData/mhr_39_20_20.txt
rm mhr