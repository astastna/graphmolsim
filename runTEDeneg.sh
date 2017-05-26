#!/bin/bash
# runs all files 001.json to 010.json for all the datasets
# $1    name of config file, has to be placed in $PWD
# $2    name of directory into which the results will be saved

# function which runs the screening
check_directories()
{
if [ ! -d $confdir/ ];then
    mkdir $confdir/
fi
if [ ! -d $confdir/$difficulty ];then
    mkdir $confdir/$difficulty
fi
if [ ! -d $confdir/$difficulty/$type ];then
    mkdir $confdir/$difficulty/$type
fi
}

run_ted_virtual_screening()
{
if [ "$i" -lt 10 ]; then
    python runTED.py -i $difficulty/random_00_30_100_20_4900/$type/00$i.json -j $difficulty/ -p $config -o $confdir/$difficulty/$type/00$i-$config-out.json
else
    python runTED.py -i $difficulty/random_00_30_100_20_4900/$type/0$i.json -j $difficulty/ -p $config -o $confdir/$difficulty/$type/0$i-$config-out.json
fi
}

# Set the needed variables according to actual location of the script
#configfile=$1
config=electronegativity
confdir=TEDeneg
#cd /storage/brno6/home/astastna/lab-export/final

# Add required modules.
module add boost-1.49
module add numpy-1.7.1-py2.7
module add scipy-0.12.0-py2.7
module add sklearn-0.14.1-py2.7

# Update environment paths.
export RDBASE=/afs/ms.mff.cuni.cz/u/s/stastna8/opt/rdkit-Release_2016_09_1
#export RDBASE=/storage/praha1/home/skodape/libs/RDKit_2016_03_2
export PYTHONPATH=$RDBASE:$PYTHONPATH
export LD_LIBRARY_PATH=$RDBASE/lib:$LD_LIBRARY_PATH


# Run the virtual screening for all the datasets
difficulty=8.0-8.5
for type in 5HT2B 5HT2C ADA2A CDK2 HDAC01 PXR_Agonist
do
    check_directories
    for i in $(seq 1 1 10)	
    do
        run_ted_virtual_screening &
    done
    wait
done

difficulty=8.5-9.0
for type in ACM1_Agonist ADA2B_Antagonist ADA2C_Antagonist CHK1
do
    check_directories
    for i in $(seq 1 1 10)	
    do
        run_ted_virtual_screening &
    done
    wait
done

difficulty=9.0-9.5
for type in 5HT1F_Agonist DRD1_Antagonist DRD2_Agonist LSHR_Antagonist OPRM_Agonist
do
    check_directories
    for i in $(seq 1 1 10)	
    do
        run_ted_virtual_screening &
    done
    wait
done

difficulty=9.8-1.0
for type in DHFR MTR1A_Agonist MTR1B_Agonist P38 V2R_Antagonist
do
    check_directories
    for i in $(seq 1 1 10)	
    do
        run_ted_virtual_screening &
    done
    wait
done

