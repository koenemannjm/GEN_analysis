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

workflowname="jeffas_dig_GEN"$kine
swif2 create $workflowname

outdirpath="/lustre19/expphy/volatile/halla/sbs/jeffas/GEN_root/simulation/GEN"$kine"/g4sbs"

preinitname=${preinit}_${nevents}events
    
for ((i=1; i<=$njobs; i++))
do

    logfile=temp/${preinitname}_dig_job_$i.txt

    outfilename=${preinitname}_job_$i.root
    outdirpathdig="/lustre19/expphy/volatile/halla/sbs/jeffas/GEN_root/simulation/GEN"$kine"/digitized"
    infilename=$outdirpath'/'$outfilename
    txtfile=${preinit}_job_$i.txt
    sbsdigjobname=${preinit}_digi_job_$i

    digscript=/w/halla-scshelf2102/sbs/jeffas/GEN_analysis/simulation/replay/run_dig.sh

    swif2 add-job -workflow $workflowname -partition production -name $sbsdigjobname -cores 1 -disk 5GB -ram 1500MB $digscript $txtfile $infilename $outdirpathdig
    #$digscript $txtfile $infilename $outdirpathdig
    #nohup $digscript $txtfile $infilename $outdirpathdig >$logfile 2>&1 &
done

#run the workflow and then print status
swif2 run $workflowname
echo -e "\n Getting workflow status.. [may take a few minutes!] \n"
swif2 status $workflowname
