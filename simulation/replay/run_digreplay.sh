SWIF_JOB_WORK_DIR=$PWD
echo 'swif_job_work_dir='$SWIF_JOB_WORK_DIR

# setup farm environments
source /site/12gev_phys/softenv.sh 2.6


# setup analyzer specific environments
export ANALYZER=/work/halla/sbs/jeffas/ANALYZER/install
source $ANALYZER/bin/setup.sh
source /work/halla/sbs/jeffas/SBS_OFFLINE/install/bin/sbsenv.sh

export SBS_REPLAY=/work/halla/sbs/jeffas/SBS_OFFLINE/SBS-replay
export ANALYZER_CONFIGPATH=$SBS_REPLAY/replay
export DB_DIR=$SBS_REPLAY/DB_MC

export OUT_DIR=$SWIF_JOB_WORK_DIR
export LOG_DIR=$SWIF_JOB_WORK_DIR

# executing the replay script
inputfile=$1
datadir=$2
outdirpath=$3

export DATA_DIR=$datadir

cp $SBS/run_replay_here/.rootrc $SWIF_JOB_WORK_DIR

analyzer -b -q 'replay_gen_mc.C+("'$inputfile'",2)'

outputfile=$OUT_DIR'/replayed_'$inputfile'.root'
mv $outputfile $outdirpath

