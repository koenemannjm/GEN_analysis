#!/bin/bash

# single-run version for swif2:

runnum=$1
maxsegments=$2

segments_per_job=1
#optional 3rd argument for output directory:
out_dir='/lustre19/expphy/volatile/halla/sbs/jeffas/GEN_root/Rootfiles'
workflowname=sbs_gen_offline

if [ $# -ge 3 ];
then
    segments_per_job=$3
fi

# very optional fourth argument: output directory:
if [ $# -ge 4 ];
then
    if [ -d $4 ]
    then
	out_dir=$4
	echo 'will write output to '$out_dir
    fi
fi

#optional workflow name
if [ $# -ge 5 ];
then
    workflowname=$5
fi


if (( $segments_per_job < 1 || $segments_per_job - 1 > $maxsegments )); then
    echo "Segments per job is not correct"
    echo $segments_per_job" Segments per job"
    exit
fi

echo 'launching swif2 replay jobs for run '$runnum', max segment = '$maxsegments', segments per job = '$segments_per_job

firstsegment=0
lastsegment=0
nsegments=0

inputstring=''

# we always start with the first segment:

# look for first segment of first stream on cache disk:
firstsegname='e1209016_'$runnum'.evio.0.0'
mssfirst='mss:/mss/halla/sbs/GEnII/raw/'$firstsegname
cachefirst='/cache/mss/halla/sbs/GEnII/raw/'$firstsegname

# no matter what we require the first stream, first segment for each job (must exist)
#testfilename='/mss/halla/sbs/GEnII/raw/'$firstsegname
# so we initialize input string to the first segment of first stream:
inputstring=' -input '$cachefirst' '$mssfirst' '

# loop over segments from zero to maxsegments
for ((i=0; i<=$maxsegments; i++))
do
    
    # the following line will be overridden by what follows:
    jobname='gen_replay_'$runnum'_segment'$i

    segment_i_added=0 # this will be set to 1 if at least segment i stream 0 is found:

    if(( $nsegments == 0 )); then # we haven't added any segments to this job yet:
	# initialize firstsegment and inputstring:
	firstsegment=$i
	inputstring=' -input '$cachefirst' '$mssfirst' '
    fi

    # look for streams zero through 2 for segment i: $i is segment number $j is stream number
    for((j=0; j<=2; j++))
    do 
	eviofilename='e1209016_'$runnum'.evio.'$j'.'$i
	mssfilename='mss:/mss/halla/sbs/GEnII/raw/'$eviofilename
	cachefile='/cache/mss/halla/sbs/GEnII/raw/'$eviofilename
    
	script='/w/halla-scshelf2102/sbs/jeffas/GEN_analysis/replay/run_GEN_swif2.sh'

	testfilename='/mss/halla/sbs/GEnII/raw/'$eviofilename
 
	# the script isn't using the following for now
	#    outfilename='match:e1209019_fullreplay_'$runnum'*seg'$i'*.root'
	#    logfilename='match:replay_gmn_'$runnum'*seg'$i'*.log'
	
	#    outcopydir=$out_dir'/rootfiles'
	#    logcopydir=$out_dir'/logs'
	
	# if the user is asking for more than one segment per job we need to concatenate all the input files together: 
	
	# behavior will depend on whether 
	
	#    addjobcm='swif2 add-job -workflow puckett_GMN_analysis -partition production'

	if [ -f "$testfilename" ]; 
	then
	    if(( $j == 0 )); then
		(( nsegments++ ))
		segment_i_added=1
	    fi
	    #echo 'Adding new swif2 job, runnum='$runnum', segment='$i 
	    if(( !($j == 0 && $i == 0) )); then # unless this is the first segment of stream 0, add the file to the job:
		inputstring+=' -input '$cachefile' '$mssfilename' '
	    fi
	fi
    done # end loop on streams 0-2
    # After adding any files found from this segment number to the job,
    # check if this is the last segment or if adding this segment caused us to reach the target number of segments per job:
    if(( $segment_i_added == 1 )); then # if we at least found stream 0 for this segment, proceed: 
	# increment the number of segments:
	# NO! This was already done above! ((nsegments++))
	if (( $i == $maxsegments || $nsegments == $segments_per_job )); 
	    then # this is either the last segment or we reached the required number of segments per job. 
	    # Go ahead and launch the job and then reset segment counter:
	    lastsegment=$i 
	    jobname='gen_replay_'$runnum'_segments'$firstsegment'_'$lastsegment
	    
	    echo 'Submitting job '$jobname' with '$nsegments' segments, runnum='$runnum
	    #		echo 'Input string = '$inputstring
	    
	    addjobcmd='swif2 add-job -workflow $workflowname -partition production -name '$jobname' -cores 1 -disk 25 GB -ram 3000 MB '$inputstring' '$script' '$runnum' -1 0 e1209016 '$firstsegment' '$nsegments' '$out_dir
	    
	    #echo 'add-job command = "'$addjobcmd'"'
	    
	    swif2 add-job -workflow $workflowname -partition production -name $jobname -cores 1 -disk 25GB -ram 3000MB $inputstring $script $runnum -1 0 e1209016 $firstsegment $nsegments $out_dir
	    
	    nsegments=0
	fi   
    else # current segment NOT added; if any segments have been added to the job, go ahead and submit the job if applicable
	if (( $nsegments > 0 ))
	then
	    ((lastsegment=$i-1))
	    jobname='gen_replay_'$runnum'_segments'$firstsegment'_'$lastsegment
	    
	    echo 'Submitting job '$jobname' with '$nsegments' segments, runnum = '$runnum
#	    echo 'Input string = "'$inputstring'"'

	    swif2 add-job -workflow $workflowname -partition production -name $jobname -cores 1 -disk 25GB -ram 3000MB $inputstring $script $runnum -1 0 e1209016 $firstsegment $nsegments $out_dir
	    
	    #addjobcmd='swif2 add-job -workflow sbs_gmn_cooking_pass1 -partition production -name '$jobname' -cores 1 -disk 25 GB -ram 3000 MB '$inputstring' '$script' '$runnum' -1 0 e1209019 '$firstsegment' '$nsegments' '$out_dir
	    #echo 'add-job command = "'$addjobcmd'"'

	    nsegments=0
	fi
	break
    fi
done
