#!/bin/bash

num_wavevectors=(10,30,100,300,500,1000)

for i in  10 30 100 300 500 600 1000 2000 20000 100000

do    
    echo "$i"
    ./eigensolver $i 1 1
done
