#!/bin/bash

num_wavevectors=(10,30,100,300,500,1000)

for i in 10 30 100 300 500 1000
do    
    echo "$i"
    ./eigensolver $i 1
    ./eigensolver $i 0
done
