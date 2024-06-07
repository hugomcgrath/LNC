#!/bin/bash

# original timestep = 10 ps
# new timestep = 20 ps
PRUNING_FACTOR=2
SYSTEMS=("LNC" "NONE")
for system in ${SYSTEMS[@]}; do
    for i in {1..20}; do
        DIR=/home/hmcgrat/LNC/$system/T$i
        XTC_IN=$DIR/traj_comp.xtc
        XTC_OUT=$DIR/traj_pruned.xtc
        echo pruning $XTC_IN
        gmx trjconv -f $XTC_IN -o $XTC_OUT -pbc atom -ur compact -skip $PRUNING_FACTOR
        echo writing pruned $XTC_OUT
    done;
done;