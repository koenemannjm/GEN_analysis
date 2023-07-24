echo 'working directory ='
echo $PWD

#SWIF_JOB_WORK_DIR=$PWD # for testing purposes
echo 'swif_job_work_dir='$SWIF_JOB_WORK_DIR

# setup farm environments
source /site/12gev_phys/softenv.sh 2.5

# Setup g4sbs specific environments
export G4SBS=/work/halla/sbs/jeffas/SBS_GEANT/g4sbs/install
source $G4SBS/bin/g4sbs.sh

# run the g4sbs command
preinit=$1
postscript=$2
nevents=$3
outfilename=$4
outdirpath=$5

rm -f $postscript

echo '/g4sbs/filename '$4 >>$postscript
echo '/g4sbs/run '$3 >>$postscript

cat $postscript

g4sbs --pre=$preinit'.mac' --post=$postscript 

# idiot proofing
if [[ ! -d $outdirpath ]]; then
    mkdir $outdirpath
fi
mv $outfilename $outdirpath
