#!/bin/sh

runs=$1;
mkdir -p ../../QC;
mkdir -p ../../Logs;
mkdir -p ../../Scripts

if [ -z "${runs}" ];
then runs=4;
fi;

rm ../../Scripts/FastQCRaw.Command*;

for i in `ls | cut -f2 -d "_"`; do echo "mkdir -p \${TMPDIR}/${i}.1 & mkdir -p \${TMPDIR}/${i}.2 & wait; lbzcat -n2 Sample_${i}/${i}_*_L*_R1_*.fastq.bz2 | fastqc -o \${TMPDIR}/${i}.1 stdin & lbzcat -n2 Sample_${i}/${i}_*_L*_R2_*.fastq.bz2 | fastqc -o \${TMPDIR}/${i}.2 stdin & wait;  cat \${TMPDIR}/${i}.1/stdin_fastqc.html | sed -e \"s:stdin:${i}_1_sequence.txt.bz2:g\" > ../../QC/${i}.1_fastqc.Raw.html & cat \${TMPDIR}/${i}.2/stdin_fastqc.html | sed -e \"s:stdin:${i}_2_sequence.txt.bz2:g\" > ../../QC/${i}.2_fastqc.Raw.html & mv \${TMPDIR}/${i}.1/stdin_fastqc.zip ../../QC/${i}.1_fastqc.Raw.zip & mv \${TMPDIR}/${i}.2/stdin_fastqc.zip ../../QC/${i}.2_fastqc.Raw.zip & wait;"; done | shuf | shuf > ../../Scripts/FastQCRaw.Command;
split ../../Scripts/FastQCRaw.Command ../../Scripts/FastQCRaw.Commands -n r/${runs};
for i in `find ../../Scripts/FastQCRaw.Commands*`; do echo "cat ${i} | parallel -j2" | qsub -l ncpus=6 -q batch -d `pwd -P` -V -N ${i}.FastQCRaw -o ../../Logs/ -e ../../Logs/; done;
exit 0;
