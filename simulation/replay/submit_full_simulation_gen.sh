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

workflowname="jeffas_sim_"$preinit
swif2 create $workflowname

outdirpath="/lustre19/expphy/volatile/halla/sbs/jeffas/GEN_root/simulation/GEN"$kine"/g4sbs"

preinitname=${preinit}_${nevents}events
    
for ((i=1; i<=$njobs; i++))
do

    outfilename=${preinitname}_job_$i.root
    postscript=${preinit}_post_$i.mac
    g4sbsjobname=${preinit}_job_$i

    g4sbsscript=/w/halla-scshelf2102/sbs/jeffas/GEN_analysis/simulation/replay/run_g4sbs.sh

    swif2 add-job -workflow $workflowname -partition production -name $g4sbsjobname -cores 1 -disk 5GB -ram 1500MB $g4sbsscript $preinit $postscript $nevents $outfilename $outdirpath
    #$g4sbsscript $preinit $postscript $nevents $outfilename $outdirpath

    outdirpathdig="/lustre19/expphy/volatile/halla/sbs/jeffas/GEN_root/simulation/GEN"$kine"/digitized"
    infilename=$outdirpath'/'$outfilename
    txtfile=${preinit}_job_$i.txt
    sbsdigjobname=${preinit}_digi_job_$i

    digscript=/w/halla-scshelf2102/sbs/jeffas/GEN_analysis/simulation/replay/run_dig.sh

    swif2 add-job -workflow $workflowname -antecedent $g4sbsjobname -partition production -name $sbsdigjobname -cores 1 -disk 5GB -ram 1500MB $digscript $txtfile $infilename $outdirpathdig
    #$digscript $txtfile $infilename $outdirpathdig

    outdirpathreplay="/lustre19/expphy/volatile/halla/sbs/jeffas/GEN_root/simulation/GEN"$kine"/replayed"
    inputfile=${preinitname}_job_$i
    digreplayjobname=${preinit}_dig_replay_job_$i

    digreplayscript=/w/halla-scshelf2102/sbs/jeffas/GEN_analysis/simulation/replay/run_digreplay.sh

    swif2 add-job -workflow $workflowname -antecedent $sbsdigjobname -partition production -name $digreplayjobname -cores 1 -disk 5GB -ram 1500MB $digreplayscript $inputfile $outdirpathdig $outdirpathreplay
    #$digreplayscript $inputfile $outdirpathdig $outdirpathreplay
done

#run the workflow and then print status
swif2 run $workflowname
echo -e "\n Getting workflow status.. [may take a few minutes!] \n"
swif2 status $workflowname
