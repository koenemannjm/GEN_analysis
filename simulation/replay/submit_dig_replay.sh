#!/bin/bash

preinit=$1
nevents=$2
njobs=$3

kine=2

# Validating the number of arguments provided
if [[ "$#" -ne 3 ]]; then
    echo -e "\n--!--\n Illegal number of arguments!!"
    echo -e " This script expects 3 arguments: <kine> <nevents> <njobs> \n"
    exit;
fi

workflowname="jeffas_dig_replay_GEN"$kine
swif2 create $workflowname

outdirpath="/lustre19/expphy/volatile/halla/sbs/jeffas/GEN_root/simulation/GEN"$kine"/g4sbs"


preinitname=${preinit}_${nevents}events
    
for ((i=21; i<=$njobs; i++))
do

    logfile=temp/${preinitname}_dig_replay_job_$i.txt

    outdirpathdig="/lustre19/expphy/volatile/halla/sbs/jeffas/GEN_root/simulation/GEN"$kine"/digitized"
    outdirpathreplay="/lustre19/expphy/volatile/halla/sbs/jeffas/GEN_root/simulation/GEN"$kine"/replayed"
    inputfile=${preinitname}_job_$i
    digreplayjobname=${preinit}_dig_replay_job_$i

    digreplayscript=/w/halla-scshelf2102/sbs/jeffas/GEN_analysis/simulation/replay/run_digreplay.sh

    swif2 add-job -workflow $workflowname -partition production -name $digreplayjobname -cores 1 -disk 5GB -ram 1500MB $digreplayscript $inputfile $outdirpathdig $outdirpathreplay
    #nohup $digreplayscript $inputfile $outdirpathdig $outdirpathreplay >$logfile 2>&1 &
    #$digreplayscript $inputfile $outdirpathdig $outdirpathreplay

done


#run the workflow and then print status
swif2 run $workflowname
echo -e "\n Getting workflow status.. [may take a few minutes!] \n"
swif2 status $workflowname
