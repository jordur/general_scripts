#!/bin/bash
procs=("$@")
# Waiting until watcher-job finishes:
while [ "`for ((i=0;i<${#procs[@]};i+=1)); do qstat | awk -v val=${procs[$i]} '{ if ($1==val) print 1 }'; done | awk '{sum=sum+$1} END {print sum}'`" > 0 ]; do
    sleep 30
done
