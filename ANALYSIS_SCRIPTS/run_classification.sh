#!/bin/bash
#===============================================================================
# Definitions for SGE
#===============================================================================
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -o /home/hmcgrat/LNC/LOGS
#$ -j y
# Find appropriate amount of memory for worst case job
#$ -l h_data=2G
# Check nuber of config files, set -t bounds to 1 to config_id_max+1
#$ -t 1-216

# Do not ignore signals:
#$ -notify
trap "" TERM

SCRIPT_PATH=$HOME/LNC/ANALYSIS_SCRIPTS
NAME=classify_selections
CONFIG_ID=$(($SGE_TASK_ID-1))
if [ $1 = "validation" ]
    then CLASSIFICATION_CONFIG_PATH=$HOME/LNC/HP_CONFIG_FILES
elif [ $1 = "test" ]
    then CLASSIFICATION_CONFIG_PATH=$HOME/LNC/HP_CONFIG_FILES_BEST
else
    echo "Invalid partition"
fi
CLASSIFICATION_CONFIG_NAME=config_$CONFIG_ID.yaml
if [ -f $CLASSIFICATION_CONFIG_PATH/$CLASSIFICATION_CONFIG_NAME ]
    then conda run -n LNC_ENV python $SCRIPT_PATH/$NAME.py $CLASSIFICATION_CONFIG_PATH $CLASSIFICATION_CONFIG_NAME
fi