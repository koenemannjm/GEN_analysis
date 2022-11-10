#!/bin/sh

Events=$1
#Jobs=$2

#declare -a Q2_range=(146 368 677 1018)
#declare -a Kine_set=("elastic" "inelastic" "beam_bkgd" "wiser_pi0" "wiser_pip" "wiser_pim")

declare -a Q2_range=(368)
declare -a Kine_set=("elastic")


for Q2 in "${Q2_range[@]}"
do

    for Kine in "${Kine_set[@]}"
    do
	
	pre_pre_script='gen/gen_'$Q2'_'$Kine'.mac'
	pre_script='job_stage/gen_'$Q2'_'$Kine'_sieve_pre.mac'
	post_script='job_stage/gen_'$Q2'_'$Kine'_sieve_post.mac'
	
	Rootfile='/lustre19/expphy/volatile/halla/sbs/jeffas/GEN_root/simulation/gen_'$Q2'_'$Kine'_'$Events'events_sieve.root'
	
	rm -f $post_script
	rm -f $pre_script
        
	cp $pre_pre_script $pre_script
	echo '/g4sbs/buildBBsieve 1' >> $pre_script	

	echo '/g4sbs/filename  '$Rootfile >> $post_script
	echo '/g4sbs/run  '$Events >> $post_script
	
	
	
	./run_gen.sh $pre_script $post_script   
	    
	    
	    

    done
done
