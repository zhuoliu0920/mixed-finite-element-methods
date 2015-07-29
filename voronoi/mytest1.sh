#!/bin/bash

k[1]=5
k[2]=10
k[3]=20
k[4]=50
k[5]=100
k[6]=200
k[7]=500
k[8]=1000
k[9]=2000
k[10]=5000

for i in 1 2 3 4 5 6 7 8 9 10 
do
echo "$i"
./a.out points="./test/uniform_dist_points_${k[i]}.txt",domain="./test/domain.txt",result="./result/uniform_dist_sin(pi*x)*sin(pi*y)_result.txt"
done
