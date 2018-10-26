#!/bin/sh

#This program takes the first name, last name and pool number for a MiSeq run and creates the file structure to run our current pipeline.  It assumes that you have already demultiplexed the whole pool, but not the individual collaborator.

umask 002

corp=$1;
collabName=$2;
pool=$3;
prefix=$4;
poolName=$5;
sampleList=$6;

if [ "${corp}" = "D" ];
then corp="DiversigenCollaborations";
else corp="CMMRCollaborations";
fi;

if [ ! -e "${sampleList}" ];
then sampleList="";
else sampleList=`readlink -e ${sampleList}`;
fi;

if [ -z "$poolName" ];
then poolName=Pool`echo $pool`;
fi;

firstName=`echo ${collabName} | cut -f1 -d "_"`;
lastName=`echo ${collabName} | cut -f2 -d "_"`;

if [ "${firstName}" = "${lastName}" ];
then firstName="";
fi;

cd /gpfs1/projects/${corp};

if [ -d "${firstName}${lastName}/${poolName}" ];
then poolName=`echo ${poolName}.Pool${pool}`;
fi;

#make the necessary directory structures
mkdir -p ${firstName}${lastName}/${poolName};
cd ${firstName}${lastName}/${poolName};
mkdir -p ${lastName}${poolName}WorkDir/Reads;
mkdir -p ${lastName}${poolName}Reads/Project_${lastName}${poolName}
mkdir -p ${lastName}${poolName}Barcodes/Project_${lastName}${poolName}/Sample_${lastName}${poolName};
mkdir -p Logs;

if [ -z "$prefix" ];
	then prefix=`echo $lastName`;
fi;

#cat in the information for files that cannot be softlinked
echo -e "FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject" > samplesheet.${lastName}${poolName}.csv
if [ -e "${sampleList}" ];
then cat /gpfs1/projects/Pools/16SV4/Pool${pool}/Pool${pool}WorkDir/SampleList | grep -wf "${sampleList}" > ${lastName}${poolName}WorkDir/SampleList;
cat /gpfs1/projects/Pools/16SV4/Pool${pool}/samplesheet.*${pool}.csv | grep -wf "${sampleList}" >> samplesheet.${lastName}${poolName}.csv
cat /gpfs1/projects/Pools/16SV4/Pool${pool}/Pool${pool}.barcodeCounts.txt | grep -wf "${sampleList}" >  ${lastName}${poolName}.barcodeCounts.txt
else cat /gpfs1/projects/Pools/16SV4/Pool${pool}/Pool${pool}WorkDir/SampleList | grep "${prefix}" > ${lastName}${poolName}WorkDir/SampleList;
cat /gpfs1/projects/Pools/16SV4/Pool${pool}/samplesheet.*${pool}.csv | grep "${prefix}" >> samplesheet.${lastName}${poolName}.csv
cat /gpfs1/projects/Pools/16SV4/Pool${pool}/Pool${pool}.barcodeCounts.txt | grep "${prefix}" >  ${lastName}${poolName}.barcodeCounts.txt
fi;

#softlink the required demultiplexed reads into the Reads/Project_* directory
for i in `ls /gpfs1/projects/Pools/16SV4/Pool${pool}/Pool${pool}Reads/Project_Pool${pool}/ | grep -f ${lastName}${poolName}WorkDir/SampleList`; do ln -sfn /gpfs1/projects/Pools/16SV4/Pool${pool}/Pool${pool}Reads/Project_Pool${pool}/$i ${lastName}${poolName}Reads/Project_${lastName}${poolName}/$i; done;

#softlink the other files in the master pool Reads directory
for i in `ls /gpfs1/projects/Pools/16SV4/Pool${pool}/Pool${pool}Reads/`; do name=`echo $i | sed "s:Pool${pool}:${lastName}${poolName}:g" | sed "s:Overall::g" | sed "s:ReagentTest::g" `; ln -sfn /gpfs1/projects/Pools/16SV4/Pool${pool}/Pool${pool}Reads/$i ${lastName}${poolName}Reads/$name; done;

#softlink the items in the individual reads into the WorkDir/Reads directory
for i in `find ${lastName}${poolName}Reads/Project_${lastName}${poolName}/Sample_*/*.bz2`; do name=`echo $i | cut -f4 -d "/" | cut -f1 -d "_"`; num=`echo $i | cut -f6 -d "_" | cut -c2`; ln -sfn ../../${i} ${lastName}${poolName}WorkDir/Reads/${name}.${num}.fq.bz2; done

#softlink the un-demultiplexed reads
for i in `ls /gpfs1/projects/Pools/16SV4/Pool${pool}/Pool${pool}Barcodes/Project_Pool${pool}/Sample_Pool${pool}/`; do name=`echo $i | sed "s:Overall::g" | sed "s:ReagentTest::g" |  sed "s:Pool${pool}:${lastName}${poolName}:g"`; ln -sfn /gpfs1/projects/Pools/16SV4/Pool${pool}/Pool${pool}Barcodes/Project_Pool${pool}/Sample_Pool${pool}/$i ${lastName}${poolName}Barcodes/Project_${lastName}${poolName}/Sample_${lastName}${poolName}/$name; done;

#softlink the files in the barcodes directory
for i in `ls /gpfs1/projects/Pools/16SV4/Pool${pool}/Pool${pool}Barcodes/ | grep -v "Project_Pool${pool}"`; do name=`echo $i | sed "s:Overall::g" | sed "s:ReagentTest::g" | sed "s:Pool${pool}:${lastName}${poolName}:g"`; ln -sfn /gpfs1/projects/Pools/16SV4/Pool${pool}/Pool${pool}Barcodes/$i ${lastName}${poolName}Barcodes/$name; done;

#softlink the log files into the Logs directory
for i in `ls /gpfs1/projects/Pools/16SV4/Pool${pool}/Logs/`; do name=`echo $i | sed "s:Overall::g" | sed "s:ReagentTest::g" |  sed "s:Pool${pool}:${lastName}${poolName}:g"`; ln -sfn /gpfs1/projects/Pools/16SV4/Pool${pool}/Logs/$i Logs/$name; done;

#softlink any remaining files into the base directory, substituting the PoolID where appropriate
for i in `ls /gpfs1/projects/Pools/16SV4/Pool${pool}/ | grep -v "Logs" | grep -v "Pool${pool}Barcodes" | grep -v "Pool${pool}Reads" | grep -v "Pool${pool}WorkDir" | grep -v "Deliverables" | grep -v "samplesheet.Pool${pool}.csv" | grep -v "barcodeCounts.txt"`; do name=`echo $i | sed "s:Overall::g" | sed "s:ReagentTest::g" |  sed "s:Pool${pool}:${lastName}${poolName}:g"`; ln -sfn /gpfs1/projects/Pools/16SV4/Pool${pool}/$i $name; done;

#echo the number of samples found, and if none are found, exit the program
numSamples=`cat ${lastName}${poolName}WorkDir/SampleList | wc -l`;
if [ $numSamples -eq 0 ];
then echo "Error: No samples found.  Please verify prefix and/or sample list.";
exit 1;
else echo "Total samples found: ${numSamples}"; 
fi;

#Allocate the threads needed for the pipeline
if [ ${numSamples} -lt 40 ];
then PROCS=$[ $[numSamples / 2] + $[numSamples % 2]];
THREADS=${numSamples};
else PROCS=20;
THREADS=40;
fi;

#automatically launch the processing job
link=`readlink -e ${lastName}${poolName}WorkDir/Reads/;`
echo "${GITREPO}/16S/fullPipelineSplit.sh ${link} ${THREADS} " | qsub -l ncpus=${PROCS} -q batch -N ${lastName}${poolName}.Process -d `pwd -P` -V;
