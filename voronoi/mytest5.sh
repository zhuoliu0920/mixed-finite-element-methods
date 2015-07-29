#!/bin/bash

k[1]=4
k[2]=16
k[3]=64
k[4]=256
k[5]=1024
k[6]=4096

for i in 3
do
echo "$i"
./a.out points="./test/uniform1_points_${k[i]}.txt",domain="./test/domain.txt",result="./result/uniform1_forFun.txt"
done
