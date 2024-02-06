#!/bin/bash

runnum=$1
first_event=$2
last_event=$3
sbs_gems=0

nevents=$(($last_event - $first_event))


fnameout_pattern='/farm_out/jeffas/gen_replayed_'$runnum'_event'$first_event'.out'
sbatch --output=$fnameout_pattern run_GEN_sbatch.csh $runnum $nevents $first_event e1209016 0 30 $sbs_gems
#run_GEN_sbatch.csh $runnum $nevents $first_event e1209016 0 30 $sbs_gems
   


