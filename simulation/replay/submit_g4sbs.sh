#!/bin/bash

preinit=$1
nevents=$2
njobs=$3

kine=2

# Validating the number of arguments provided
if [[ "$#" -ne 3 ]]; then
    echo -e "\n--!--\n Illegal number of arguments!!"
    echo -e " This script expects 3 arguments: <preinit> <nevents> <njobs> \n"
    exit;
fi

workflowname="jeffas_g4sbs_"$preinit
swif2 create $workflowname

outdirpath="/lustre19/expphy/volatile/halla/sbs/jeffas/GEN_root/simulation/GEN"$kine"/g4sbs"

preinitname=${preinit}_${nevents}events
    
for ((i=1; i<=$njobs; i++))
do

    logfile=temp/${preinitname}_job_$i.txt

    outfilename=${preinitname}_job_$i.root
    postscript=${preinit}_post_$i.mac
    g4sbsjobname=${preinit}_job_$i

    g4sbsscript=/w/halla-scshelf2102/sbs/jeffas/GEN_analysis/simulation/replay/run_g4sbs.sh

    #swif2 add-job -workflow $workflowname -partition production -name $g4sbsjobname -cores 1 -disk 5GB -ram 1500MB $g4sbsscript $preinit $postscript $nevents $outfilename $outdirpath
    $g4sbsscript $preinit $postscript $nevents $outfilename $outdirpath
    #nohup $g4sbsscript $preinit $postscript $nevents $outfilename $outdirpath >$logfile 2>&1 &
done

#run the workflow and then print status
swif2 run $workflowname
echo -e "\n Getting workflow status.. [may take a few minutes!] \n"
swif2 status $workflowname
