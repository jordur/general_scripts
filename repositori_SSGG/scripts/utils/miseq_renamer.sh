#!/bin/bash
if [ $# != 3 ] ; then
    echo '###################################'
    echo '#### MiSeq sample renamer 0.1 #####'
    echo '###################################'
    echo 'Argument 1: MiSeq run directory'
    echo 'Arugment 2: Initial number value'
    echo 'Argument 3: Number value you want to use for replacement'
    echo 'This will rename all sample .fastq.gz files inside every sample directory from 001 to specified number.'
    echo 'ON1-17281_S6_L001_R1_001.fastq.gz -> ON1-17281_S6_L001_R1_00X.fastq.gz'
    echo ''
    echo 'Usage: miseq_renamer.sh /path/to/miseq/runfolder number_to_rename'
    echo 'ie: miseq_renamer.sh /share/gluster/ToDo/MiSeq/130731_M01167_0022_000000000-A5G7K 1 3'
    echo ''
    exit 0
fi

NameOfRun=$1
InitialCount=$2
NewCount=$3


cd $NameOfRun
for directory in $(find * -type d ) ; do
        cd $directory
#               for i in $(find -iname "*_00${InitialCount}.fastq.gz" | sed 's/\.\/\(.\+\)_00${InitialCount}\.fastq\.gz/\1/') ; do
                for i in $(find -iname "*_00${InitialCount}.fastq.gz" | sed "s/\.\/\(.\+\)_00${InitialCount}\.fastq\.gz/\1/") ; do
                        mv ${i}_00${InitialCount}.fastq.gz ${i}_00${NewCount}.fastq.gz
        done
        cd ..
done
