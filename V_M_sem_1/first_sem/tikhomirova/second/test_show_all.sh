#!/bin/bash
for((i=0;i<101;i++))
do
./a.out 30 $i
done |grep "answer:"|less
