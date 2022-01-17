#!/bin/bash
# vim:foldmethod=marker
# a script that predicts outcomes of single patients, re-trains 
# the model periodically and assesses the accuracy of the predictions

# SETUP {{{
WORK_DIR="/home/cluster/mlabar/data/realtime"

SEQLIB="${WORK_DIR}/realtime_seqlib_2.fas"                          # sequences are kept in this FASTA-file
CC_SCRIPT="/home/cluster/mlabar/bin/cluster-compare2.py"         # Python-script that generates the dataset(s)
R_SCRIPT="/home/cluster/mlabar/bin/realtime_modeling.R"          # R-script that generates the predictive model

HIVTRACE_T="0.01"   # HIV-Trace distance threshold for clustering
HIVTRACE_MO="500"   # HIV-Trace minimum overlap of sequences
HIVTRACE_FR="0.05"  # HIV-Trace maximum ambiguity fraction to resolve

INTERV=3        # (integer) Number of years over which the outcome is calculated
USE_N=5000      # (integer) the last N sequences to be used in a training dataset
LAB_INTERV=60	# (integer) number of days that are allowed between sequence timestamp and sample_day
#}}}

# INPUT {{{
while [ -n "$1" ]; do
    case "$1" in
        -n) NEW_FASTA="$(realpath $2)"                          # a FASTA-file with the patient's sequence
            NEW_HEADER="$(grep "^>" $NEW_FASTA)"
            NEW_PID="$(echo $NEW_HEADER | cut -f 2 -d "|")"     # patient_id
            NEW_TIME="$(echo $NEW_HEADER | cut -f 3 -d "|")"    # timestamp
            #NEW_PATIENT_DATA="$(realpath $3)"                  # a CSV-file with predictors (UNUSED)

            echo "--- (realtime.sh) Timestamp: $NEW_TIME ---"
            echo "--- (realtime.sh) Patient ID: $NEW_PID ---"
            ;;

        -t) OUTCOME_INT="$2"
            ;;

    esac
    shift
done
#}}}

# MAIN LOOP {{{
cd $WORK_DIR

# append new sequence
cat $NEW_FASTA >> $SEQLIB

# only continue if the patient is known
if [ $((NEW_PID)) -eq 999999 ]; then
    echo "--- (realtime.sh) Patient ID is 999999. Exiting. ---"
    exit
fi

# calculate past timestamps
if [[ "$(hostname)" =~ "marco-mbp" ]]; then
    P2="$(date -j -f "%FT%H:%M:%S" "$NEW_TIME" "+%F")"          # sequence timestamp (only date)
else
    P2="$(date -d "$NEW_TIME" "+%F")"                           # sequence timestamp (only date)
fi

# cluster with HIV-Trace
echo "--- (realtime.sh) Clustering with HIV-Trace ---"
hivtrace -i $SEQLIB -a resolve -r HXB2_pol -t $HIVTRACE_T -m $HIVTRACE_MO -g $HIVTRACE_FR

# process results into tables
echo "--- (realtime.sh) Analysing clusters ---"
$CC_SCRIPT -f $SEQLIB -j ${SEQLIB}.results.json -d $P2 -p "$NEW_PID"

# train models and make predictions
echo "--- (realtime.sh) Model fitting and prediction ---"
Rscript "$R_SCRIPT" "$P2" "$NEW_PID" "$USE_N" "$LAB_INTERV"
#}}}
