#!/bin/bash

k[1]=2
k[2]=8
k[3]=32
k[4]=128
k[5]=512
k[6]=2048

for i in 1 2 3 4 5 6
do
echo "$i"
./a.out points="./test/uniform3_points_${k[i]}.txt",domain="./test/domain.txt",result="./result/uniform3_sin(pi*x)*sin(pi*y)_result.txt"
done
