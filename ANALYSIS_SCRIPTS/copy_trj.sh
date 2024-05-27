#!/bin/bash

ORIGIN_DIR=/netmount/ribfit1/ribfit1/ribosome_antibiotics/compare_w_wo_antibiotics/analyse/rmsd
DESTINATION_DIR=/home/hmcgrat/LNC
rsync -uvh --progress $ORIGIN_DIR/pdbs/r50S_lnc_minimized_ren.pdb $DESTINATION_DIR/LNC/system.pdb
rsync -uvh --progress $ORIGIN_DIR/pdbs/r50S_wo_antibiotic_minimized_ren.pdb $DESTINATION_DIR/NONE/system.pdb
for i in {1..20}; do
    echo LNC $i
    rsync -uvh --progress $ORIGIN_DIR/xtcs/r50S_lnc_${i}_equil_20-120ns_compact.xtc $DESTINATION_DIR/LNC/T${i}/traj_comp.xtc
    echo NONE $i
    rsync -uvh --progress $ORIGIN_DIR/xtcs/r50S_wo_antibiotic_${i}_equil_20-120ns_compact.xtc $DESTINATION_DIR/NONE/T${i}/traj_comp.xtc
done