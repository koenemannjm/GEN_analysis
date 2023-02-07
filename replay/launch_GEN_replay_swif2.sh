#!/bin/bash

# single-run version for swif2:

kinematic=$1
runnum=$2
firstseg=$3
maxsegments=$4

segments_per_job=1

if [ $# -eq 4 ];
then
    segments_per_job=$4
fi

echo 'launching swif2 replay jobs for run '$runnum', max segment = '$maxsegments', segments per job = '$segments_per_job

for ((i=$firstseg; i<=$maxsegments; i++))
do
    fnameout_pattern='/farm_out/jeffas/gen_replayed_'$runnum'_segment'$i'.out'
    #    sbatch --output=$fnameout_pattern run_GMN_sbatch_nohodo.sh $runnum -1 0 e1209019 $i 1
    jobname='gen_replay_'$runnum'_segment'$i
    
    # look for first segment on cache disk:
    firstsegname='e1209016_'$runnum'.evio.0.0'
    mssfirst='mss:/mss/halla/sbs/GEnII/raw/'$firstsegname
    cachefirst='/cache/mss/halla/sbs/GEnII/raw/'$firstsegname
    
    eviofilename1='e1209016_'$runnum'.evio.0.'$i
    mssfilename1='mss:/mss/halla/sbs/GEnII/raw/'$eviofilename1
    cachefile1='/cache/mss/halla/sbs/GEnII/raw/'$eviofilename1

    eviofilename2='e1209016_'$runnum'.evio.1.'$i
    mssfilename2='mss:/mss/halla/sbs/GEnII/raw/'$eviofilename2
    cachefile2='/cache/mss/halla/sbs/GEnII/raw/'$eviofilename2

    eviofilename3='e1209016_'$runnum'.evio.2.'$i
    mssfilename3='mss:/mss/halla/sbs/GEnII/raw/'$eviofilename3
    cachefile3='/cache/mss/halla/sbs/GEnII/raw/'$eviofilename3
    
    script='/work/halla/sbs/jeffas/GEN_analysis/replay/run_GEN_swif2.sh'

    testfilename1='/mss/halla/sbs/GEnII/raw/'$eviofilename1
    testfilename3='/mss/halla/sbs/GEnII/raw/'$eviofilename3
 
    # the script isn't using the following for now
    #    outfilename='match:e1209019_fullreplay_'$runnum'*seg'$i'*.root'
    #    logfilename='match:replay_gmn_'$runnum'*seg'$i'*.log'
    
    #    outcopydir=$out_dir'/rootfiles'
    #    logcopydir=$out_dir'/logs'

    # if the user is asking for more than one segment per job we need to concatenate all the input files together: 
    

    if [ -f "$testfilename3" ]; 
    then
	echo 'Adding new swif2 job, runnum='$runnum', segment='$i 
    
	if [ $i -gt 0 ]
	then
	    echo 'segment '$i' also requires first segment'
	    swif2 add-job -workflow jeffas_${kinematic}_analysis -partition production -name $jobname -cores 1 -disk 25GB -ram 4000MB -input $cachefile1 $mssfilename1 -input $cachefile2 $mssfilename2 -input $cachefile3 $mssfilename3 -input $cachefirst $mssfirst $script $runnum -1 0 e1209016 $i 1
	else
	    echo 'segment '$i' IS first segment'
	    swif2 add-job -workflow jeffas_${kinematic}_analysis -partition production -name $jobname -cores 1 -disk 25GB -ram 4000MB -input $cachefile1 $mssfilename1 -input $cachefile2 $mssfilename2 -input $cachefile3 $mssfilename3 -input $cachefirst $mssfirst $script $runnum -1 0 e1209016 $i 1
	fi
    elif [ -f "$testfilename1" ]; 
    then
	echo 'Adding new swif2 job, runnum='$runnum', segment='$i
	
	if [ $i -gt 0 ]
	then
	    echo 'segment '$i' also requires first segment'
	    swif2 add-job -workflow jeffas_${kinematic}_analysis -partition production -name $jobname -cores 1 -disk 25GB -ram 4000MB -input $cachefile1 $mssfilename1 -input $cachefirst $mssfirst $script $runnum -1 0 e1209016 $i 1
	else
	    echo 'segment '$i' IS first segment'
	    swif2 add-job -workflow jeffas_${kinematic}_analysis -partition production -name $jobname -cores 1 -disk 25GB -ram 4000MB -input $cachefile1 $mssfilename1 -input $cachefirst $mssfirst $script $runnum -1 0 e1209016 $i 1
	fi
    fi
done
