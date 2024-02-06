#!/bin/bash

#SBATCH --partition=production
#SBATCH --account=halla
#SBATCH --mem-per-cpu=1500

#cd /work/halla/sbs/puckett/GMN_ANALYSIS

echo 'working directory ='
echo $PWD

#SWIF_JOB_WORK_DIR=$PWD
echo 'swif_job_work_dir='$SWIF_JOB_WORK_DIR



MODULES=/etc/profile.d/modules.sh 
if [[ $(type -t module) != function && -r ${MODULES} ]]; then 
    source ${MODULES} 
fi 
module use /group/halla/modulefiles
module load analyzer/1.7.4
module list

# setup environment for ANALYZER and SBS-offline:

echo 'working directory = '$PWD

export SBS=/work/halla/sbs/jeffas/SBS_OFFLINE/install
source $SBS/bin/sbsenv.sh

export SBS_REPLAY=/work/halla/sbs/jeffas/SBS_OFFLINE/SBS-replay
export DB_DIR=$SBS_REPLAY/DB
export DATA_DIR=/cache/mss/halla/sbs/GEnII/raw

export OUT_DIR=$SWIF_JOB_WORK_DIR
export LOG_DIR=$SWIF_JOB_WORK_DIR
#export OUT_DIR="/volatile/halla/sbs/jeffas/GEN_root/Rootfiles"
#export LOG_DIR="/volatile/halla/sbs/jeffas/GEN_root/logs"

echo 'OUT_DIR='$OUT_DIR
echo 'LOG_DIR='$LOG_DIR

#export OUT_DIR=/volatile/halla/sbs/puckett/GMN_REPLAYS/SBS4/rootfiles
#export LOG_DIR=/volatile/halla/sbs/puckett/GMN_REPLAYS/SBS4/rootfiles
export ANALYZER_CONFIGPATH=$SBS_REPLAY/replay

runnum=$1
maxevents=$2
firstevent=$3

prefix=$4
firstsegment=$5
maxsegments=$6

outcopydir=$7

cp /work/halla/sbs/jeffas/GEN_analysis/replay/.rootrc $SWIF_JOB_WORK_DIR

#currently set up for no SBS GEM analysis and no CM plots
analyzer -b -q 'replay_gen.C+('$runnum','$maxevents','$firstevent','\"$prefix\"','$firstsegment','$maxsegments',2,0,0,0)'

outfilename=$OUT_DIR'/e1209016_*'$runnum'*.root'
logfilename=$LOG_DIR'/replay_gen_'$runnum'*.log' 

mv $outfilename $outcopydir


