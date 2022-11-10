#!/bin/sh

declare -a Q2_range=(29 66 97)
declare -a Sigmas=(3 2.5 2)



for Q2 in "${Q2_range[@]}"
do

    for sigma in "${Sigmas[@]}"
    do

	root -b -q 'ftrglog_trgrt_gmn.C('$Q2','$sigma')'
    #root -b -q 'hc_tgrt_beam_gmn.C('$Q2')'

    done
done
