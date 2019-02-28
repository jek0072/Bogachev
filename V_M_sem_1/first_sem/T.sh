#!/bin/bash
for ((n=2; n<=150;n++)) ; do for ((m=3;m<=n;m+=3)) ; do echo "n=$n m=$m" ; ./a.out $n $m 1; done ; done 
