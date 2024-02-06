#!/bin/bash

runnum=$1
first_event=$2
last_event=$3
sbs_gems=0

nevents=$(($last_event - $first_event))

workflowname=jeffas_sbs_offline
outdir=/lustre19/expphy/volatile/halla/sbs/jeffas/GEN_root/Rootfiles/GEM_luminosity/zero_suppressed
jobname=gen_replay_${runnum}_event${first_event}

swif2 create $workflowname

swif2 add-job -workflow $workflowname -partition production -name $jobname -cores 1 -disk 25GB -ram 3000MB /w/halla-scshelf2102/sbs/jeffas/GEN_analysis/replay/run_GEN_swif2.sh $runnum $nevents $first_event e1209016 0 30 $outdir

   
swif2 run $workflowname
