#!/bin/bash

export inputFiles=$1;
#Input is a list of directories that need to be merged.
export outputDirectory=$2;
export THREADS=$3;
if [ -z "${THREADS}" ];
	then export THREADS=`grep -c ^processor /proc/cpuinfo`;
fi


#Check to see if the merge directory already exists
if [ ! -d $outputDirectory ];
then mkdir $outputDirectory;
echo "Made output directory at $outputDirectory";
fi

export outputDirectory=`readlink -e \`echo $outputDirectory\``;
export outname=`echo $outputDirectory | rev | cut -f1 -d "/" | rev`
echo "Output directory is now $outputDirectory";

#Makes the directory structure for the merged directory
mkdir -p ${outputDirectory}/Logs;
mkdir -p ${outputDirectory}/${outname}Reads/Project_${outname};
mkdir -p ${outputDirectory}/${outname}WorkDir/Reads;
mkdir -p ${outputDirectory}/${outname}Barcodes/Project_${outname}/Sample_${outname};
echo "FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject" > ${outputDirectory}/samplesheet.${outname}.csv

#Cycles thru each input directory and links/merges the needed read files for each sample
for i in `cat ${inputFiles}`;
do export input=`readlink -e ${i}`;
echo "Adding pool ${input} to the merged set"
#Build the samplesheet.csv and check for barcode collisions
for j in `cat ${input}/samplesheet.*.csv | grep -f ${input}/*WorkDir/SampleList`; do barcode=`echo ${j} | cut -f5 -d ","`;
echo ${barcode};
name=`echo ${j} | cut -f3 -d ","`;
if [ -z ${TMPDIR}/BarcodeList ];
then echo ${barcode} > ${TMPDIR}/BarcodeList;
elif cat ${TMPDIR}/BarcodeList | grep -q "${barcode}";
then testname=`cat ${outputDirectory}/samplesheet.${outname}.csv | grep ${barcode} | cut -f3 -d ","`;
if [ "${name}" != "${testname}" ];
then echo -e "Colliding barcode found: ${barcode}\n Affected samples:\n${name}\n${testname}\n";
rm -rf ${outputDirectory} && exit 1;
fi;
else
echo ${j} >> ${outputDirectory}/samplesheet.${outname}.csv;
echo ${barcode} >> ${TMPDIR}/BarcodeList;
fi;
done;
#Create the files for the merged reads
find ${input}/*Reads/Project_* | grep -v ".bz2" | grep -v ".csv" | grep "Sample_" | grep -f ${input}/*WorkDir/SampleList | parallel -j${THREADS} -I {} 'name=`echo "{}" | rev | cut -f1 -d "/" | rev`;
if [ ! -d ${TMPDIR}/${outname}Reads/Project_${outname}/${name} ];
then mkdir -p ${TMPDIR}/${outname}Reads/Project_${outname}/${name};
echo "Made directory ${TMPDIR}/${outname}Reads/Project_${outname}/${name}";
else echo "Merging files in directory ${TMPDIR}/${outname}Reads/Project_${outname}/${name}";
fi;
for file in `ls {} | grep -v ".csv"`;
do fileName=`echo ${file} | sed "s:.bz2::g" | sed -e "s:_[ACGNT]\{12\}_:_NNNNNNNNNNNN_:g"`;
echo "Adding ${file} to the merged file ${fileName}";
bzcat {}/${file} >> ${TMPDIR}/${outname}Reads/Project_${outname}/${name}/${fileName};
done;
for file in `ls {} | grep ".csv"`;
do cat {}/${file} >> ${TMPDIR}/${outname}Reads/Project_${outname}/${name}/${file};
done;';
#Make the Monolithic Fastq Files
find ${input}/*Barcodes/Project_*/Sample_* | grep ".gz" | parallel -j${THREADS} -I {} 'name=`echo {} | rev | cut -f1 -d "/" |rev | cut -f2- -d "_" | sed "s:.gz::g"`;
if [ ! -d ${TMPDIR}/${outname}Barcodes/Project_${outname}/Sample_${outname} ];
then mkdir -p ${TMPDIR}/${outname}Barcodes/Project_${outname}/Sample_${outname};
echo "Made directory ${TMPDIR}/${outname}Barcodes/Project_${outname}/Sample_${outname}";
else echo "Merging monolithic fastqs in directory ${TMPDIR}/${outname}Barcodes/Project_${outname}/Sample_${outname}";
fi;
echo "Merging the undemultiplexed reads and barcodes file {} into ${outname}Barcodes/Project_${outname}/Sample_${outname}/${outname}_${name}";
zcat {} >> ${TMPDIR}/${outname}Barcodes/Project_${outname}/Sample_${outname}/${outname}_${name};' &
#Cat together the required stats files to get accurate read statistics
cat ${input}/sampleSheet.notDemultiplexed.*.csv >> ${outputDirectory}/sampleSheet.notDemultiplexed.${outname}.csv &
wait;
done;

#Zip up the required files after all have been integrated
echo "Compressing demultiplexed fastq files for use in the pipeline";
find ${TMPDIR}/${outname}Reads/Project_${outname}/Sample_*/*.fastq | parallel -j${THREADS} -I {} '
bzip2 {}' &
wait;
cp -r ${TMPDIR}/${outname}Reads/Project_${outname}/* ${outputDirectory}/${outname}Reads/Project_${outname}/;

cd ${TMPDIR}/${outname}Barcodes/Project_${outname}/Sample_${outname}/;
echo "Compressing the undemultiplexed fastq files for use in the pipeline";
ls | grep ".fastq" | parallel -j${THREADS} -I {} '
pigz {} && cp -r {}.gz ${outputDirectory}/${outname}Barcodes/Project_${outname}/Sample_${outname}/' &
wait;

#Creates simlinks to the read files in the working directory
cd ${outputDirectory}/${outname}WorkDir/Reads;
for j in `ls ${outputDirectory}/${outname}Reads/Project_${outname}`
do name=`echo $j | rev | cut -f1 -d "/" | rev | sed 's:Sample_::g'`;
for seq in `ls ${outputDirectory}/${outname}Reads/Project_${outname}/${j} | grep "bz2"`
do link=`echo $seq | sed -r 's:_[ACGNT]+_L001_R1_001.fastq:.1.fq:g' | sed -r 's:_[ACGNT]+_L001_R2_001.fastq:.2.fq:g'`;
ln -s ${outputDirectory}/${outname}Reads/Project_${outname}/${j}/${seq} ${outputDirectory}/${outname}WorkDir/Reads/${link};
done;
done;
#Make the new SampleList
ls  ${outputDirectory}/${outname}WorkDir/Reads/ | grep ".1.fq.bz2" | sed "s:.1.fq.bz2::g" | sort | uniq > ../SampleList
cd ${outputDirectory};
for i in `cat ${outputDirectory}/${outname}WorkDir/SampleList`; do count=`bzcat ${outputDirectory}/${outname}WorkDir/Reads/${i}.1.fq.bz2 | wc -l`; count=`echo $((count / 4))`; echo -e "${i}\t${count}" >> ${outputDirectory}/${outname}.barcodeCounts.txt;
done;
rm -rf ${TMPDIR}/*;
${GITREPO}/16S/fullPipelineSplit.sh ${outputDirectory}/${outname}WorkDir/Reads/ ${THREADS};
exit 0;
