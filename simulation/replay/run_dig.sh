#!/bin/bash

# setup farm environments
source /site/12gev_phys/softenv.sh 2.6

# Setup g4sbs specific environments
export LIBSBSDIG=/work/halla/sbs/jeffas/SBS_GEANT/libsbsdig/install
source $LIBSBSDIG/bin/sbsdigenv.sh

dbfile=$LIBSBSDIG/db/db_gmn_conf.dat


# run the sbsdig command
txtfile=$1 # .txt file containing input file paths
infilename=$2
outdirpath=$3

echo $infilename >>$txtfile


sbsdig $dbfile $txtfile

mv $infilename $outdirpath
