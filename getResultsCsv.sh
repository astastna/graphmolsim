#!/bin/sh
# This program gets the results from folders for all datasets
# and prints a summary csv file to stdout.
#
# Parameters:
# $1    algorithm from which to get the data - mcs, ged, ted, s-mcs, s-ged
# $2    config file - mcs:  config-Elements-Any.json
#                   - ged:  Beam-0.8 (file ends with .prop, but don't give it into the parameter)
# $3    parent folder of the given run (where all the data are saved)
# $4    type of output file - mcs: out / murcko (in output file typically: 
#                                               001-config-Elements-Any.json-murcko.json)
#                           - ged: GED / murcko
# $5    what should be evaluated: auc, efXX where XX are hundredths: 005, 01, 02, 05

print_header(){
    # Prints the first line of the csv table
    for i in 01 02 03 04 05 06 07 08 09 10
    do
        printf ";0%s" "$i"
    done
    printf "; avg \n"
}

results_to_csv(){
    # Processes one dataset -> prints one line of resulting csv
    # Goes through the directories, gets the results form the json files
    # and prints them in form of a csv table.
    # Last column of the table is average of all results for the given dataset.
    sum=0
    avg=0
    printf "%s;" "$dataset"
    # For all splits:
    for i in 01 02 03 04 05 06 07 08 09 10
    do
        # Get the result from the correct output file
        case $algtype in
            mcs)
                content="./$dir/$difficulty/$dataset/0$i-$configfile-$outtype.json"
                ;;
            ged)
        	alpha=$( cat $content | grep '"'alpha'"' | sed 's/^.*://' | sed 's/,.*$//')
        	match=$( cat $content | grep '"'matching'"' | sed 's/^.*://' | sed 's/,.*$//')
                content="./$dir/$difficulty/$dataset/0$i-$match-$alpha-$outtype.json"
                ;;
            ted)
                ;;
            s-mcs)
                ;;
            s-ged)
                ;;
            * ) 
                echo "Unknown algorithm type:" $1
                echo "Possible algorithms are: mcs, ged, ted, s-mcs, s-ged."
                exit 1 
                ;;
        esac

        case $whattoeval in
            auc) 
                value=auc
                ;;
            ef005)
                value=0.005
                ;;
            ef01)
                value=0.01
                ;;
            ef02)
                value=0.02
                ;;
            ef05)
                value=0.05
                ;;
            *)
                echo "Unknown evaluation type:" $1
                echo "Possible things to evaluate are: auc, ef005, ef01, ef02, ef05."
                exit 1
                ;;
        esac
 
        # Take line with value and 
        # remove all but the result number
        result=$( cat $content | grep '"'$value'"' | sed 's/^.*://' | sed 's/,.*$//')
        # Print the result with cell separator at the end to prepare the next cell
        printf "%s;" "$result"
        # Add the result to the sum of results for this dataset
        sum=`echo $sum + $result | bc -l`
    done
    # Count the average result
    avg=`echo $sum / 10.0 | bc -l`
    printf "0%s \n" "$avg"
}

# Parsing the arguments from command line

if [ -n "${1}" ]; then
    algtype=$1
else
    echo "Parameter 1: Type of algorithm has to be set."
    exit 1
fi

if [ -n "${2}" ]; then
    configfile=$2
else
    echo "Parameter 2: Config file has to be given."
    exit 1
fi

if [ -n "${3}" ]; then
    dir=$3
else
    echo "Parameter 3: Directory with program data has to be given."
    exit 1
fi

if [ -n "${4}" ]; then
    whattoeval=$4
else
    echo "Parameter 5: Choose what should be evaluated: auc, ef005, ef01, ef02, ef05."
    exit 1
fi

# Go through all datasets and process them,
# with difficulty change print a line with difficulty
print_header
difficulty=8.0-8.5
echo $difficulty
for dataset in 5HT2B 5HT2C ADA2A CDK2 HDAC01 PXR_Agonist
do
    results_to_csv
done

difficulty=8.5-9.0
echo $difficulty 
for dataset in ACM1_Agonist ADA2B_Antagonist ADA2C_Antagonist CHK1 
do
    results_to_csv
done

difficulty=9.0-9.5
echo $difficulty 
for dataset in 5HT1F_Agonist DRD1_Antagonist DRD2_Agonist LSHR_Antagonist OPRM_Agonist 
do
    results_to_csv
done

difficulty=9.8-1.0
echo $difficulty 
for dataset in DHFR MTR1A_Agonist MTR1B_Agonist P38 V2R_Antagonist 
do
    results_to_csv
done
