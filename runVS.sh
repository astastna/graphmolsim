#!/bin/sh
# This script should be run in the root folder where are 
# - datasets in folders according to difficulty
# - directory config with the configuration files
# - files convertToGXL.py, convertToCXL.py, virtual-screening-configurable.py, graph-edit-distance.jar
#
# This script runs the virtual screening according to the given parameters.
# Individual splits are run in paralell if possible. (Max. 10 CPU's used.)
#
# Parameters:
# $1    algorithm from which to get the data - mcs, ged, ted, s-mcs, s-ged
# $2    config file - mcs:  config-Elements-Any.json
#                   - ged:  config-eneg.json
# $3    parent folder of the given run (where all the data are saved)
# $4    type of output file - mcs: out / murcko (in output file typically: 
#                                               001-config-Elements-Any.json-murcko.json)
#                           - ged: GED / murcko

###################################################################
# Functions start here                                            #
###################################################################

create_gxl(){
    # The number of split doesn't matter, because there are same files in all
    # config files for the splits.
    datanum=001 
    
python convertToGXL.py -i $difficulty/random_00_30_100_20_4900/$type/$datanum.json -j $difficulty/ -o GXLfiles/$difficulty/$type/
}
###################################################################

check_directories()
{
if [ ! -d $dir/ ];then
    mkdir $dir/
fi
if [ ! -d $dir/$difficulty ];then
	mkdir $dir/$difficulty
fi
if [ ! -d $dir/$difficulty/$type ];then
	mkdir $dir/$difficulty/$type/
fi
if [ ! -d GXLfiles/ ];then
    mkdir GXLfiles/
fi
if [ ! -d GXLfiles/$difficulty ];then
	mkdir GXLfiles/$difficulty
fi
if [ ! -d GXLfiles/$difficulty/$type ];then
	mkdir GXLfiles/$difficulty/$type/
fi
}
###################################################################

run_ged_virtual_screening()
{
# Set the number of dataset (we want 00i but not 0010)
if [ "$i" -lt 10 ]; then
    datanum=00$i
else
    datanum=0$i
fi

# Check that rest of the directories exist, create them if not
if [ ! -d $dir/$difficulty/$type/$datanum/ ]; then
	mkdir $dir/$difficulty/$type/$datanum/
fi
if [ ! -d GXLfiles/$difficulty/$type/molecules/ ]; then
	mkdir GXLfiles/$difficulty/$type/molecules/
fi

# Create cxl files for the GED call
python convertToCXL.py -i $difficulty/random_00_30_100_20_4900/$type/$datanum.json -j $difficulty/ -o $dir/$difficulty/$type/$datanum

# Prepare data to create config file for the call
alg=$(cat config/$configfile | grep matching | sed 's/^.*: *"//' | sed 's/".*$//') # This doesn't have to be taken from the configuration file 
alpha=$(cat config/$configfile | grep alpha | sed 's/^.*: *"//' | sed 's/".*$//')  # In order to get results for more combinations easily, just make a for loop through alg and alpha 
									  	   # function ged
sourcePath="source=$datapath/$dir/$difficulty/$type/$datanum/source.cxl"
targetPath="target=$datapath/$dir/$difficulty/$type/$datanum/target.cxl"
path="path=$datapath/GXLfiles/$difficulty/$type/molecules/"
output="result=$datapath/$dir/$difficulty/$type/$datanum/"
matching="matching="$alg                      # algorithm used for the matching       
node="node="$(cat config/$configfile | grep nodeInsertDel | sed 's/^.*: *//' | sed 's/,.*$//')
edge="edge="$(cat config/$configfile | grep edgeInsertDel | sed 's/^.*: *//' | sed 's/,.*$//')
node0importance="nodeAttr0Importance="$(cat config/$configfile | grep nodeValence | sed 's/^.*: *//' | sed 's/ *,.*$//')
node1importance="nodeAttr1Importance="$(cat config/$configfile | grep nodeElement | sed 's/^.*: *//' | sed 's/ *,.*$//')
node2importance="nodeAttr2Importance="$(cat config/$configfile | grep nodeElectronegativity | sed 's/^.*: *//' | sed 's/ *,.*$//')
node3importance="nodeAttr3Importance="$(cat config/$configfile | grep nodePharmacophores | sed 's/^.*: *//' | sed 's/ *,.*$//')
edge0importance="edgeAttr0Importance="$(cat config/$configfile | grep edgeValence | sed 's/^.*: *//' | sed 's/ *,.*$//')
propfile=$(cat config/$configfile | grep fullConfiguration | sed 's/^.*: *"//' | sed 's/".*$//')
# fullConfiguration has to be the last one in the list of configurations,
# otherwise there will be a bug here


# Create the new config file by substituing the original values
# node - cost of node deletion / insert  
# edge - cost of edge deletion / insert  
# node0attribute - valence                         
# node1attribute - element                         
# node2attribute - electronegativity               
# node3attribute - pharmacophores
# edge0attribute - valence                         

awk -v SOURCE=$sourcePath -v TARGET=$targetPath -v GXLPATH=$path -v RESULT=$output -v MATCHING=$matching -v NODE=$node -v EDGE=$edge -v NODE0IMPORTANCE=$node0importance -v NODE1IMPORTANCE=$node1importance -v NODE2IMPORTANCE=$node2importance -v NODE3IMPORTANCE=$node3importance -v EDGE0IMPORTANCE=$edge0importance -v ALPHA="alpha="$alpha '{
    sub(/^source=.*$/, SOURCE);                
    sub(/^target=.*$/, TARGET);
    sub(/^path=.*$/, GXLPATH);
    sub(/^result=.*$/, RESULT);
    sub(/^matching=.*$/, MATCHING);
    sub(/^node=.*$/, NODE);
    sub(/^edge=.*$/, EDGE);
    sub(/^nodeAttr0Importance=.*$/, NODE0IMPORTANCE);
    sub(/^nodeAttr1Importance=.*$/, NODE1IMPORTANCE);
    sub(/^nodeAttr2Importance=.*$/, NODE2IMPORTANCE);
    sub(/^nodeAttr3Importance=.*$/, NODE3IMPORTANCE);
    sub(/^edgeAttr0Importance=.*$/, EDGE0IMPORTANCE);
    sub(/^alpha=.*$/, ALPHA);
    print;
}' config/$propfile > $dir/$difficulty/$type/$datanum/$alg-$alpha-$datanum.prop

# Run the GED counting
java -jar graphEditDistance.jar $dir/$difficulty/$type/$datanum/$alg-$alpha-$datanum.prop

# Get the results and create the file with AUC and data for further evaluation
python evaluate-GED-screening.py -i $difficulty/random_00_30_100_20_4900/$type/$datanum.json -j $difficulty/ -o $dir/$difficulty/$type/$datanum-$alg-$alpha-$outtype.json -g $dir/$difficulty/$type/$datanum/$alg-$alpha.json
}
###################################################################

run_mcs_virtual_screening()
{
if [ "$i" -lt 10 ]; then
    python virtual-screening-configurable.py -i $difficulty/random_00_30_100_20_4900/$type/00$i.json -j $difficulty/ -o $dir/$difficulty/$type/00$i-$configfile-out.json -c config/$configfile
else
    python virtual-screening-configurable.py -i $difficulty/random_00_30_100_20_4900/$type/0$i.json -j $difficulty/ -o $dir/$difficulty/$type/0$i-$configfile-out.json -c config/$configfile
fi
}
###################################################################

ged(){
    # Check, that all needed directories exist.
    # This is here as there are no split-number dependent directories.
    check_directories
    
    # Create gxl files for the given dataset (time-consuming, so we won't do it 10-times)
    # Removed because we are going to use GXLfiles folder with precomputed gxl files.
    # Can be changed in run_ged_virtual_screening() when generating property file 
    # (GXLfiles -> $dir)
     create_gxl
     for i in $(seq 1 1 10)
          do
          run_ged_virtual_screening &
          done
     wait # wait here is important, because $type, $difficulty are global
      # and so we can't run two datasets at once
      # can be repaired by giving these variables as arguments of the function
      # but when running with 10 CPUs it is not an issue
}
###################################################################

mcs(){
    check_directories
    for i in $(seq 1 1 10)
        do
            run_mcs_virtual_screening &
        done
    wait
}
###################################################################

run_dataset(){
    # Run the dataset on correct algorithm
        case $algtype in
            mcs)
                mcs
                ;;
            ged)
                ged
                ;;
            ted)
                ;;
            s-mcs)
                ;;
            s-ged)
                ;;
            * )
                echo "Unknown algorithm type:" $algtype
                echo "Possible algorithms are: mcs, ged, ted, s-mcs, s-ged."
                exit 1
                ;;
        esac
}
###################################################################
# The main program starts here                                    #
###################################################################

# Update environment paths.
#export RDBASE=/afs/ms.mff.cuni.cz/u/s/stastna8/opt/rdkit-Release_2016_09_1
#export RDBASE=/storage/praha1/home/skodape/libs/RDKit_2016_03_2
export PYTHONPATH=$RDBASE:$PYTHONPATH
export LD_LIBRARY_PATH=$RDBASE/lib:$LD_LIBRARY_PATH

# Add required modules.
#module add boost-1.49
#module add numpy-1.7.1-py2.7
#module add scipy-0.12.0-py2.7
#module add sklearn-0.14.1-py2.7
#module add jdk-8

# Save the type of run
algtype=$1
configpath=$2
configfile=$(echo $configpath | sed 's#.*/##') # remove the part before last /
dir=$3 # TODO make a parameter
outtype=out
datapath=.

# Run the virtual screening for all the datasets
difficulty=8.0-8.5
for type in 5HT2C ADA2A CDK2 HDAC01 PXR_Agonist 5HT2B
do
   run_dataset
done
wait

difficulty=8.5-9.0
for type in ACM1_Agonist ADA2B_Antagonist ADA2C_Antagonist CHK1
do
   run_dataset
done


difficulty=9.0-9.5
for type in 5HT1F_Agonist DRD1_Antagonist DRD2_Agonist LSHR_Antagonist OPRM_Agonist
do
   run_dataset
done

difficulty=9.8-1.0
for type in DHFR MTR1A_Agonist MTR1B_Agonist P38 V2R_Antagonist
do
   run_dataset
done
