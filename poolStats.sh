#!/bin/sh

THREADS=$1;

if [ -z "$THREADS" ];
        then THREADS=`grep -c ^processor /proc/cpuinfo`;
fi

DIR=`pwd -P | cut -f7 -d "/" | cut -f1 | sed 's/WorkDir//g'`;
totalReads=`cat Logs/${DIR}.phixOverall.stats.txt | grep "reads" | cut -f1 -d " ";`
phiXReads=`cat Logs/${DIR}.phixOverall.stats.txt | grep "exactly" | cut -f5 -d " "`;
phiXReads1=`cat Logs/${DIR}.phixOverall.stats.txt | grep ">1" | cut -f5 -d " "`;
phiXReads=$((${phiXReads}+${phiXReads1}));
barcoded=0;
for i in `cat ${DIR}.barcodeCounts.txt | cut -f2`; do barcoded=$(($barcoded + `echo $i`)); done;
strictMerged=`cat Deliverables/${DIR}.Stats.StrictMerge.Combined.txt | grep "Merged Stats" -A2 | grep "Total" | cut -f2 -d " "`;
standardMerged=`cat Deliverables/${DIR}.Stats.StandardMerge.Combined.txt | grep "Merged Stats" -A2 | grep "Total" | cut -f2 -d " "`;
numSamples=`cat ${DIR}.barcodeCounts.txt | wc -l`;
totalReads=$((${totalReads} * 250));
barcoded=$((${barcoded} * 500));
phiXReads=$((${phiXReads} * 250));
strictMerged=$((${strictMerged} * 500));
standardMerged=$((${standardMerged} * 500));
percentBarcoded=$((100*${barcoded}/${totalReads}));
percentPhiX=$((100*${phiXReads}/${totalReads}));
percentStrict=$((100*${strictMerged}/${barcoded}));
percentStandard=$((100*${standardMerged}/${barcoded}));
echo -e "PoolID\t#Samples\tMachine Reads\tBarcoded Reads\t%Barcoded Reads:Machine Reads\t#phiX Mapped\t%phiX:Machine Reads\tStrict Merged Reads\t%Strict Merged Reads:Barcoded Reads\tStandard Merged Reads\t%Standard Merged Reads:Barcoded Reads\n${DIR}\t${numSamples}\t${totalReads}\t${barcoded}\t${percentBarcoded}%\t${phiXReads}\t${percentPhiX}%\t${strictMerged}\t${percentStrict}%\t${standardMerged}\t${percentStandard}%" > Deliverables/${DIR}.PoolStats.txt;
