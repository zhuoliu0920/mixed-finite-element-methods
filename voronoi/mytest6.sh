#!/bin/bash

k[1]=4
k[2]=16
k[3]=64
k[4]=256
k[5]=1024
k[6]=4096

for i in 1 2 3 4 5 6
do
echo "$i"
./a.out points="./test/uniform2_points_${k[i]}.txt",domain="./test/domain.txt",result="./result/uniform2_sin(pi*x)*sin(pi*y)_result.txt"
done
