#!/bin/bash
#
# This is the main script.
# It runs the script which performes the virtual screening and then 
# runs another script to get csv's with summary of the results.
#
# This script should be run in the root folder where are 
# - datasets in folders according to difficulty
# - directory config with the configuration files
# - files convertToGXL.py, convertToCXL.py, virtual-screening-configurable.py, graph-edit-distance.jar
#
# Parameters:
# $1    config file - mcs:  config-Elements-Any.json
#                   - ged:  config-eneg.json
# [$2]    parent folder of the given run (where all the data are saved)

check_params(){
if [ -n "${1}" ]
then
    config=$1
    configfile=$(echo $config | sed 's#.*/##') # remove the part before last /

else
	echo " Not enough arguments given. Please pass following parameters to the script:"
	echo " * path to the config file (choose one in config folder or create your own)"
	echo " * (optional) name of output folder (default one containing name of config file will be used if not given)"
	exit 0
fi
}

set_dir(){
if [ -n "${1}" ] ; then
	dir=$1
else
	dir=$alg-$configfile
fi
}

show_algorithm_choice(){
echo "Choose graph algorithm:"
echo "m) Maximum common (edge) subgraph"
echo "g) Graph edit distance"
echo "or exit (e)."
}

read_options(){
	local choice
	read -p "Type m, g or e and pres enter: " choice
	case $choice in
		m) alg=mcs
		   ;;
		g) alg=ged
		   ;;
		e) exit 0;;
		*) echo -e "${RED} Wrong key pressed ${STD}" && sleep 2
	esac
}

# Check and save the parameters
check_params $1

# Show dialog
show_algorithm_choice
read_options

# Set output directory
set_dir $2

# Run the virtual screening:
./runVS.sh $alg $config $dir

# Evaluate the results:
for what in auc ef005 ef01 ef02 ef05
do
./getResultsCsv.sh  $alg $confdir $dir > results/$alg-$confdir-$what.csv
done
