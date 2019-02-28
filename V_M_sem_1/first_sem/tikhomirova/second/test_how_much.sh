#!/bin/bash
for((i=0;i<101;i++))
do
./a.out 100 $i
done |grep "answer:"|wc
