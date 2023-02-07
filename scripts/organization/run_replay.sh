#!/bin/bash

kinematic=$1
runnum=$2


cd ../replay

launch_GEN_replay_swif2.sh $kinematic $runnum 0 30
#echo "launch_GEN_replay_swif2.sh "$kinematic" "$runnum" 0 30"

