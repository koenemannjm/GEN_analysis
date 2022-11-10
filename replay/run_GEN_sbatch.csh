#!/bin/tcsh

#SBATCH --partition=production
#SBATCH --account=halla
#SBATCH --chdir=/work/halla/sbs/puckett/GMN_ANALYSIS
#SBATCH --mem-per-cpu=4000
#SBATCH --constraint=farm19

source /site/12gev_phys/softenv.csh 2.5
cd /work/halla/sbs/jeffas/GEN_analysis/replay


# setup environment for ANALYZER and SBS-offline:

#setenv ANALYZER /work/halla/sbs/ANALYZER/install
setenv ANALYZER /work/halla/sbs/jeffas/ANALYZER/install
source $ANALYZER/bin/setup.csh
#source /work/halla/sbs/SBS_OFFLINE/install/bin/sbsenv.csh
#source /work/halla/sbs/jeffas/SBS_OFFLINE/install/bin/sbsenv.csh
source /work/halla/sbs/jeffas/SBS_OFFLINE/install/bin/sbsenv.csh

setenv SBS_REPLAY /work/halla/sbs/jeffas/SBS_OFFLINE/SBS-replay
#setenv SBS_REPLAY /work/halla/sbs/jeffas/SBS-replay
setenv DB_DIR $SBS_REPLAY/DB
setenv DATA_DIR /cache/halla/sbs/raw
#setenv DATA_DIR /volatile/halla/sbs/jeffas/GMN_root/raw
setenv OUT_DIR /volatile/halla/sbs/jeffas/GEN_root/Rootfiles
setenv LOG_DIR /volatile/halla/sbs/jeffas/GEN_root/logs
setenv ANALYZER_CONFIGPATH $SBS_REPLAY/replay

set runnum=$1
set maxevents=$2
set firstevent=$3

set prefix=$4
set firstsegment=$5
set maxsegments=$6



#set cmd="analyzer "${ANALYZER_CONFIGPATH}"'/replay_gmn_nohodo.C+("$runnum","$maxevents","$firstevent","\"$prefix\"","$firstsegment","$maxsegments")'"

#echo $cmd

#analyzer -q
analyzer  -q 'replay_gen.C+('$runnum','$maxevents','$firstevent','\"$prefix\"','$firstsegment','$maxsegments',2,0,1,0)'
#analyzer  -q 'replay_gen_helicity_test.C+('$runnum','$maxevents','$firstevent','\"$prefix\"','$firstsegment','$maxsegments',2,0,0,0,0)'
#analyzer  -q 'replay_all_GEMs.C+('$runnum','$firstsegment','$maxsegments','\"$prefix\"','$firstevent','$maxevents',1)'
#analyzer  -q 'replay_SBSGEM.C+('$runnum','$maxevents','$firstevent','\"$prefix\"','$firstsegment','$maxsegments')'



