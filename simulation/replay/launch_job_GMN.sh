#!/bin/sh

Events=$1
Jobs=$2

#declare -a Q2_range=(146 368 677 1018)
#declare -a Kine_set=("elastic" "inelastic" "beam_bkgd" "wiser_pi0" "wiser_pip" "wiser_pim")

declare -a Q2_range=(136)
declare -a Kine_set=("elastic" "inelastic" "wiser_pi0" "wiser_pip" "wiser_pim")


for Q2 in "${Q2_range[@]}"
do

    for Kine in "${Kine_set[@]}"
    do
	
	pre_script='gmn/gmn_'$Q2'_'$Kine'.mac'
	
	for((i=1; i<=$Jobs; i++))
	do
	    
	    pre_pre_script='gmn/gmn_'$Q2'_'$Kine'.mac'
            pre_script='job_stage/gmn_'$Q2'_'$Kine'_pre_job'$i'.mac'
            post_script='job_stage/gmn_'$Q2'_'$Kine'_post_job'$i'.mac'
	    
	    Rootfile='/lustre19/expphy/volatile/halla/sbs/jeffas/GMN_root/simulation/gmn_'$Q2'_'$Kine'_sieve_'$Events'events_job'$i'.root'
	   
            rm -f $post_script
            rm -f $pre_script
	    
            cp $pre_pre_script $pre_script
            echo '/g4sbs/buildBBsieve 1' >> $pre_script
	    
            echo '/g4sbs/filename  '$Rootfile >> $post_script
            echo '/g4sbs/run  '$Events >> $post_script
	    
	    
	    rm -f $post_script
            
	    echo '/g4sbs/filename  '$Rootfile >> $post_script
	    echo '/g4sbs/run  '$Events >> $post_script
	    
	    
	    fnameout_pattern='/farm_out/jeffas/gmn_simulation_'$Q2'_'$Kine'_job'$i'.out'
	    sbatch --output=$fnameout_pattern ./run_gmn.sh $pre_script $post_script   
	    #./run_gmn.sh $pre_script $post_script   
	    
	    
	    
	done
    done
done
