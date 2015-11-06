#!/bin/sh

#This program takes the first name, last name and pool number for a MiSeq run and creates the file structure to run our current pipeline.  It assumes that you have already demultiplexed the whole pool, but not the individual collaborator.

firstName=$1;
lastName=$2;
pool=$3;

mkdir ${firstName}${lastName};
mkdir ${firstName}${lastName}/Pool${pool};
cd ${firstName}${lastName}/Pool${pool};
mkdir ${lastName}Pool${pool}WorkDir;
mkdir ${lastName}Pool${pool}WorkDir/Reads;
cat ../../StatsProject/16S/Pool${pool}/Pool${pool}WorkDir/SampleList | grep "${lastName}" > ${lastName}Pool${pool}WorkDir/SampleList;
for i in `ls ../../StatsProject/16S/Pool${pool}/Pool${pool}WorkDir/Reads/ | grep -f ${lastName}Pool${pool}WorkDir/SampleList | grep "bz2"`; do ln -s ../../../../StatsProject/16S/Pool${pool}/Pool${pool}WorkDir/Reads/$i ${lastName}Pool${pool}WorkDir/Reads/$i; done;
mkdir ${lastName}Pool${pool}Reads;
mkdir ${lastName}Pool${pool}Reads/Project_${lastName}Pool${pool}
for i in `ls ../../StatsProject/16S/Pool${pool}/Pool${pool}Reads/Project_Pool${pool}/ | grep -f ${lastName}Pool${pool}WorkDir/SampleList`; do ln -s ../../../../StatsProject/16S/Pool${pool}/Pool${pool}Reads/Project_Pool${pool}/$i ${lastName}Pool${pool}Reads/Project_${lastName}Pool${pool}/$i; done;
for i in `ls ../../StatsProject/16S/Pool${pool}/Pool${pool}Reads/`; do name=`echo $i | sed "s:Pool${pool}:${lastName}Pool${pool}:g" | sed "s:Overall::g" | sed "s:ReagentTest::g" `; ln -s ../../../StatsProject/16S/Pool${pool}/Pool${pool}Reads/$i ${lastName}Pool${pool}Reads/$name; done;
mkdir ${lastName}Pool${pool}Barcodes;
mkdir ${lastName}Pool${pool}Barcodes/Project_${lastName}Pool${pool};
mkdir ${lastName}Pool${pool}Barcodes/Project_${lastName}Pool${pool}/Sample_${lastName}Pool${pool};
for i in `ls ../../StatsProject/16S/Pool${pool}/Pool${pool}Barcodes/Project_Pool${pool}/Sample_Pool${pool}/`; do name=`echo $i | sed "s:Overall::g" | sed "s:ReagentTest::g" |  sed "s:Pool${pool}:${lastName}Pool${pool}:g"`; ln -s ../../../../../StatsProject/16S/Pool${pool}/Pool${pool}Barcodes/Project_Pool${pool}/Sample_Pool${pool}/$i ${lastName}Pool${pool}Barcodes/Project_${lastName}Pool${pool}/Sample_${lastName}Pool${pool}/$name; done;
for i in `ls ../../StatsProject/16S/Pool${pool}/Pool${pool}Barcodes/ | grep -v "Project_Pool${pool}"`; do name=`echo $i | sed "s:Overall::g" | sed "s:ReagentTest::g" | sed "s:Pool${pool}:${lastName}Pool${pool}:g"`; ln -s ../../../StatsProject/16S/Pool${pool}/Pool${pool}Barcodes/$i ${lastName}Pool${pool}Barcodes/$name; done;
mkdir Logs;
for i in `ls ../../StatsProject/16S/Pool${pool}/Logs/`; do name=`echo $i | sed "s:Overall::g" | sed "s:ReagentTest::g" |  sed "s:Pool${pool}:${lastName}Pool${pool}:g"`; ln -s ../../../StatsProject/16S/Pool${pool}/Logs/$i Logs/$name; done;
for i in `ls ../../StatsProject/16S/Pool${pool}/ | grep -v "Logs" | grep -v " Pool${pool}Barcodes" | grep -v "Pool${pool}Reads" | grep -v "Pool${pool}WorkDir" | grep -v "Deliverables"`; do name=`echo $i | sed "s:Overall::g" | sed "s:ReagentTest::g" |  sed "s:Pool${pool}:${lastName}Pool${pool}:g"`; ln -s ../../StatsProject/16S/Pool${pool}/$i $name; done;
cd ${lastName}Pool${pool}WorkDir/Reads/;
echo "~gesell/Programs/ToTest/fullPipelineSplit.sh" | qsub -l ncpus=20 -q batch -N ${lastName}Pool${pool}.Process -d `pwd -P` -V;
