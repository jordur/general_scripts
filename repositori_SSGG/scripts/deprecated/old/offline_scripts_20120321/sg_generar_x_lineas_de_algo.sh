#!/bin/bash
yoquiero=$1
lineas=$2
let sol=$lineas+1
myvar=1
until [ $myvar -eq $sol ]
do
    echo $1
    myvar=$(( $myvar + 1 ))
done
