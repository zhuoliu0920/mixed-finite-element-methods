#!/bin/bash

k[1]=4
k[2]=16
k[3]=64

n[1]=3
n[2]=5
n[3]=9

p[1]=0.5
p[2]=1
p[3]=2
p[4]=5
p[5]=7.5
p[6]=10

for i in 1 2 3 4 5 6
do
	for j in 1 2 3
	do
		echo "-----------------------------------------------------------------"
		echo "${k[j]} points and ${n[j]} layers and ${p[i]} perturbation"
		./0constant.out points="./test/uniform1_points_${k[j]}.txt",domain="./test/domain.txt",result="./result_perturb_new/0constant_result.txt",nlayer=${n[j]},percent=${p[i]}
		./1linear.out points="./test/uniform1_points_${k[j]}.txt",domain="./test/domain.txt",result="./result_perturb_new/1linear_result.txt",nlayer=${n[j]},percent=${p[i]}
		./2square.out points="./test/uniform1_points_${k[j]}.txt",domain="./test/domain.txt",result="./result_perturb_new/2square_result.txt",nlayer=${n[j]},percent=${p[i]}
		./3product.out points="./test/uniform1_points_${k[j]}.txt",domain="./test/domain.txt",result="./result_perturb_new/3product_result.txt",nlayer=${n[j]},percent=${p[i]}
		./4cubic.out points="./test/uniform1_points_${k[j]}.txt",domain="./test/domain.txt",result="./result_perturb_new/4cubic_result.txt",nlayer=${n[j]},percent=${p[i]}
		./5sincos.out points="./test/uniform1_points_${k[j]}.txt",domain="./test/domain.txt",result="./result_perturb_new/5sincos_result.txt",nlayer=${n[j]},percent=${p[i]}
	done
done
