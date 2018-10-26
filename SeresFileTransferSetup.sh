#!/bin/sh

THREADS=$1;
link=$2;

if [ -z "${THREADS}" ];
then export THREADS=`grep -c ^processor /proc/cpuinfo`;
fi

if [ -z "${link}" ];
then link=`pwd -P`;
else link=`readlink -e ${link}`;
fi;

Project=`echo ${link} | cut -f5,6 -d "/" | tr -s "/" "."`;
directory=`echo ${Project} | tr -s "." "_" | sed -e "s:Batch:Batch_:g"`;
if [[ ${Project} == *".Pool"* ]];
then isRerun=1;
Project=`echo ${Project} | sed -e "s:Pool[0-9]*:RR:g"`;
directory=`echo ${directory}| sed -e "s:Pool[0-9]*:RR:g"`;
fi;
batch=`echo ${Project} | cut -f2 -d "."`;

mkdir -p ${link}/Deliverables/${directory};
ln -s ${link}/Deliverables/Raw_Read1.fq.bz2 ${link}/Deliverables/${directory}/${batch}_read1.fastq.bz2;
ln -s ${link}/Deliverables/Raw_Read2_Barcodes.fq.bz2 ${link}/Deliverables/${directory}/${batch}_barcodes.fastq.bz2;
ln -s ${link}/Deliverables/Raw_Read3.fq.bz2 ${link}/Deliverables/${directory}/${batch}_read2.fastq.bz2;
cat ${link}/Deliverables/Demultiplex_Sheet.txt | sed -e "s:Seres\.::g" > ${link}/Deliverables/${directory}/Demultiplex_Sheet.txt
echo -e "Client\tVendor\tRefence Agreement Version\tTransfer Date\tNames of Files Transferred\tFile Size" > ${link}/Deliverables/${directory}/${Project}.Manifest.txt;
for i in `ls ${link}/Deliverables/${directory} | grep -v "${Project}.MD5Sums.txt\|${Project}.Manifest.txt"`; do echo "sum=\`md5sum ${link}/Deliverables/${directory}/${i} | cut -f1 -d ' '\`; size=\`du -shL ${link}/Deliverables/${directory}/${i} | cut -f1\`; echo -e \"Seres Therapudic\tDiversigen\tWMS DTA Agreement for study ID SER-262-001 dated 28_APR_2017\t\`date '+%D'\`\t${i}\t\${size}\" >> ${link}/Deliverables/${directory}/${Project}.Manifest.txt & echo -e \"\${sum}\t${i}\" >> ${link}/Deliverables/${directory}/${Project}.MD5Sums.txt;"; done | parallel -j${THREADS} -k;
exit 0;
