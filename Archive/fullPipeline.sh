#!/bin/sh

THREADS=$1;


if [ -z "$THREADS" ];
	then THREADS=`grep -c ^processor /proc/cpuinfo`;
fi
DIR=`pwd -P | cut -f9 -d "/" | cut -f1 | sed 's/WorkDir//g'`;
echo ${DIR};
cat ../SampleList | parallel -j${THREADS} -I {} "echo {}";
echo "Temporary Directory: ${TMPDIR}";
cat ../SampleList | parallel -j${THREADS} -I {} "bzcat {}.1.fq.bz2 | sed 's/1:N:0:.*/1:N:0:/g' > ${TMPDIR}/{}.1.fq";
cat ../SampleList | parallel -j${THREADS} -I {} "bzcat {}.2.fq.bz2 | sed 's/2:N:0:.*/3:N:0:/g' > ${TMPDIR}/{}.2.fq";
cat ../SampleList | parallel -j${THREADS} -I {} "usearch70 -fastq_mergepairs ${TMPDIR}/{}.1.fq -reverse ${TMPDIR}/{}.2.fq -fastq_minovlen 50 -fastq_maxdiffs 4 -fastq_truncqual 5 -fastqout ${TMPDIR}/{}.MergedRaw.fq";
cat ../SampleList | parallel -j${THREADS} -I {} "echo {}; usearch70 -fastq_filter ${TMPDIR}/{}.MergedRaw.fq -fastq_maxee .5 -fastqout ${TMPDIR}/{}.filteredRaw.fq";
cat ../SampleList | parallel -j${THREADS} -I {} "usearch70 -fastq_mergepairs ${TMPDIR}/{}.1.fq -reverse ${TMPDIR}/{}.2.fq -fastq_minovlen 50 -fastq_maxdiffs 0 -fastq_minmergelen 252 -fastq_maxmergelen 254 -fastq_truncqual 5 -fastqout ${TMPDIR}/{}.merged.fq";
cat ../SampleList | parallel -j${THREADS} -I {} "usearch70 -fastq_filter ${TMPDIR}/{}.merged.fq -fastq_maxee .05 -relabel "${TMPDIR}/{}_" -fastqout ${TMPDIR}/{}.filtered.fq";
cat ${TMPDIR}/*.filtered.fq > ${TMPDIR}/seqs.fq.temp;
cat ${TMPDIR}/*.filteredRaw.fq > ${TMPDIR}/seqs.raw.fq.temp;
cat ${TMPDIR}/seqs.fq.temp | sed 's:@/space1/tmp/[0-9]*.sug-moab/:@:g' > ${TMPDIR}/seqs.fq;
cat ${TMPDIR}/seqs.raw.fq.temp | sed 's:@/space1/tmp/[0-9]*.sug-moab/:@:g' > ${TMPDIR}/seqs.raw.fq;
bowtie2 -x /users/mcwong/work/references/phix/phix -U ${TMPDIR}/seqs.fq --end-to-end --very-sensitive --reorder -p ${THREADS} --un ${TMPDIR}/seqs.filtered.fq -S /dev/null 2>../../Logs/phix.bleed.txt;
bowtie2 -x /users/mcwong/work/references/phix/phix -U ${TMPDIR}/seqs.raw.fq --end-to-end --very-sensitive --reorder -p ${THREADS} --un ${TMPDIR}/seqs.raw.filtered.fq -S /dev/null 2>../../Logs/phix.raw.prep.bleed.txt;
mkdir ../split_libraries;
fq2fa ${TMPDIR}/seqs.filtered.fq ${TMPDIR}/seqs.fna;
mv ${TMPDIR}/seqs.fna ../split_libraries;
mv ${TMPDIR}/seqs.raw.filtered.fq ${TMPDIR}/MergedRaw.fq;
cat ../SampleList | parallel -j1 -I {} "cat ${TMPDIR}/{}.1.fq" > ${TMPDIR}/Read1.fq;
cat ../SampleList | parallel -j1 -I {} "cat ${TMPDIR}/{}.2.fq" > ${TMPDIR}/Read2.fq;
fq2fa ${TMPDIR}/MergedRaw.fq ${TMPDIR}/MergedRaw.fa;
cat ${TMPDIR}/Read1.fq ${TMPDIR}/Read2.fq > ${TMPDIR}/Reads.fq;
fq2fa ${TMPDIR}/Reads.fq ${TMPDIR}/Reads.fa;
usearch70 -fastq_mergepairs ${TMPDIR}/Read1.fq -reverse ${TMPDIR}/Read2.fq -fastaout ${TMPDIR}/temp.fa;
cat ${TMPDIR}/temp.fa >> ${TMPDIR}/Reads.fa
usearch70 -usearch_global ${TMPDIR}/MergedRaw.fa -db ~mcwong/jobs/Francisella/Francisella_V4.udb -strand plus -id .968 -uc ${TMPDIR}/FrancisellaMerged.uc -maxaccepts 0 -maxrejects 0 -threads ${THREADS};
usearch70 -usearch_global ${TMPDIR}/Reads.fa -db ~mcwong/jobs/Francisella/Francisella_V4.udb -strand both -id .968 -uc ${TMPDIR}/FrancisellaRaw.uc -maxaccepts 0 -maxrejects 0 -threads ${THREADS};
cat ${TMPDIR}/FrancisellaMerged.uc | cut -f9,10 | grep -v "*$" | cut -f1 | cut -f1 -d " " > ${TMPDIR}/Remove;
cat ${TMPDIR}/FrancisellaRaw.uc | cut -f9,10 | grep -v "*$" | cut -f1 | cut -f1 -d " " >> ${TMPDIR}/Remove;
cat ${TMPDIR}/MergedRaw.fq | perl /users/mcwong/fastqfilter.pl -v ${TMPDIR}/Remove > ${TMPDIR}/Temp.fq;
mv ${TMPDIR}/Temp.fq ${TMPDIR}/MergedRaw.fq;
cat ${TMPDIR}/Read1.fq | perl /users/mcwong/fastqfilter.pl -v ${TMPDIR}/Remove > ${TMPDIR}/Temp.fq;
mv ${TMPDIR}/Temp.fq ${TMPDIR}/Read1.fq;
cat ${TMPDIR}/Read2.fq | perl /users/mcwong/fastqfilter.pl -v ${TMPDIR}/Remove > ${TMPDIR}/Temp.fq;
mv ${TMPDIR}/Temp.fq ${TMPDIR}/Read2.fq;
pbzip2 -c ${TMPDIR}/Read1.fq > ${TMPDIR}/Read1.fq.bz2;
pbzip2 -c ${TMPDIR}/Read2.fq > ${TMPDIR}/Read2.fq.bz2;
pbzip2 -c ${TMPDIR}/MergedRaw.fq > ${TMPDIR}/MergedRaw.fq.bz2;
cat ../SampleList | xargs -I {} mkdir {};
cat ../SampleList | parallel -j${THREADS} -I {} 'cat ${TMPDIR}/{}.2.fq | perl /users/mcwong/fastqfilter.pl -v ${TMPDIR}/Remove | pbzip2 -p1 -c > {}/{}.2.fq.bz2';
cat ../SampleList | parallel -j${THREADS} -I {} 'cat ${TMPDIR}/{}.1.fq | perl /users/mcwong/fastqfilter.pl -v ${TMPDIR}/Remove | pbzip2 -p1 -c > {}/{}.1.fq.bz2';
mkdir ${TMPDIR}/uparse;
usearch70 -derep_fulllength ../split_libraries/seqs.fna  -output ${TMPDIR}/uparse/derep.fna -sizeout -uc ${TMPDIR}/uparse/derep.uc;
usearch70 -sortbysize       ${TMPDIR}/uparse/derep.fna -output ${TMPDIR}/uparse/sorted.fa -minsize 2;
cp ${TMPDIR}/uparse/sorted.fa ${TMPDIR}/uparse/temp.fa
for i in {0.4,0.8,1.2,1.6,2.0,2.4,2.8,3.2};
do 
usearch70 -cluster_otus     ${TMPDIR}/uparse/temp.fa -otus   ${TMPDIR}/uparse/temp1.fa -otu_radius_pct $i -uc ${TMPDIR}/uparse/cluster_$i.uc -fastaout ${TMPDIR}/uparse/clustering.$i.fasta.out;
cat ${TMPDIR}/uparse/clustering.$i.fasta.out | grep "^>" | grep chimera | sed 's/^>//g' | sed -re 's/;n=.*up=/\t/g' | sed 's/;$//g' | tee -a ${TMPDIR}/uparse/chimeras.txt > ${TMPDIR}/uparse/chimeras.$i.txt;
cat ${TMPDIR}/uparse/clustering.$i.fasta.out | grep "^>" > ${TMPDIR}/uparse/uparseref.decisions.$i.txt;
rm ${TMPDIR}/uparse/clustering.$i.fasta.out;
mv ${TMPDIR}/uparse/temp1.fa ${TMPDIR}/uparse/temp.fa;
done
mv ${TMPDIR}/uparse/temp.fa ${TMPDIR}/uparse/otus1.fa
usearch70 -uchime_ref       ${TMPDIR}/uparse/otus1.fa  -db /users/mcwong/work/references/gold/gold.fa -strand plus -nonchimeras ${TMPDIR}/uparse/otus.fa -uchimeout ${TMPDIR}/uparse/uchimeref.uc;
usearch70 -usearch_global   ${TMPDIR}/uparse/otus.fa -db /stornext/snfs7/cmmr/dpsmith/silva/latest/silva_V4.udb -id .968 -strand plus -threads $THREADS -uc ${TMPDIR}/uparse/otus2taxa.uc -maxaccepts 0 -maxrejects 0;
cat ${TMPDIR}/uparse/derep.fna | grep -A1 "size=1;" | cut -f2 -d ">" | getSeq ${TMPDIR}/uparse/derep.fna > ${TMPDIR}/uparse/singletons.fna;
usearch70 -usearch_global ${TMPDIR}/uparse/singletons.fna -db ${TMPDIR}/uparse/sorted.fa -id .99 -uc ${TMPDIR}/uparse/singletons2otus.uc -strand plus -threads $THREADS -maxaccepts 32 -maxrejects 128 -minqt 1 -leftjust -rightjust -wordlength 12;
resolveIterativeUparse.pl ${TMPDIR}/uparse/cluster_*.uc ${TMPDIR}/uparse/singletons2otus.uc ${TMPDIR}/uparse/otus2taxa.uc --derep ${TMPDIR}/uparse/derep.uc --chimeras ${TMPDIR}/uparse/chimeras.txt --uchime ${TMPDIR}/uparse/uchimeref.uc --taxonomy /stornext/snfs7/cmmr/dpsmith/silva/latest/silva.map;
mv otu_table.biom ${TMPDIR}/uparse/;
biom summarize-table -i ${TMPDIR}/uparse/otu_table.biom -o ${TMPDIR}/uparse/stats.otu_table.txt;
cat ../split_libraries/seqs.fna | grep "^>" | cut -f1 -d "_" | cut -f2 -d ">" | sort | uniq -c > ../Stats.MergedReads.txt;
cat ${TMPDIR}/uparse/stats.otu_table.txt | tail -n +17 | sed 's/^ //g' | sed -re 's/: /\t/g' | sed 's/\.0$//g' > ../Stats.MappedReads.txt;
perl ~mcwong/gitRepo/16S_workflows/StatsComparisonMergedVsMapped.pl ../Stats.MergedReads.txt ../Stats.MappedReads.txt > ../Stats.Combined.txt
tar -cvjf ../uparse.tar.bz2 ${TMPDIR}/uparse;
~mcwong/gitRepo/16S_workflows/recoverBarcodesForRaw.pl ${TMPDIR}/Read1.fq.bz2 ../../${DIR}Barcodes/Project_${DIR}/Sample_${DIR}/${DIR}_NoIndex_L001_R2_001.fastq.gz | pbzip2 -c > ${TMPDIR}/RawReadsBarcodes.fq.bz2;
~mcwong/gitRepo/16S_workflows/recoverBarcodesForRaw.pl ${TMPDIR}/MergedRaw.fq.bz2 ../../${DIR}Barcodes/Project_${DIR}/Sample_${DIR}/${DIR}_NoIndex_L001_R2_001.fastq.gz | pbzip2 -c > ${TMPDIR}/MergedBarcodes.fq.bz2;
mkdir ../../Deliverables;
for i in `ls ${TMPDIR}/*.fq.bz2`; do name=`echo $i | cut -f5 -d "/"`; cp $i ../../Deliverables/${DIR}.${name};done
summarize_taxa.py -i ${TMPDIR}/uparse/otu_table.biom -o ../../Deliverables/${DIR}.TaxaSummary;
cp ${TMPDIR}/uparse/otu_table.biom ../../Deliverables/${DIR}.otu_table.biom; 
cp ../Stats.Combined.txt ../../Deliverables/${DIR}.Stats.Combined.txt;
cp ../SampleList ../../Deliverables/${DIR}.SampleList;
cat ../../samplesheet.${DIR}.csv | cut -f3,5 -d "," | tr "," "\t" > ../../Deliverables/${DIR}.SampleSheet.txt
head -1 ~mcwong/IlluminaHeaderExample > ../../Deliverables/${DIR}.ExampleQiimeMappingFile.txt
tail -n+2 ../../Deliverables/${DIR}.SampleSheet.txt | sed -re 's/(.*)\t(.*)/\1\t\2\tGGACTACHVGGGTWTCTAAT\tGTGCCAGCMGCCGCGGTAA\t\1/g' >> ../../Deliverables/${DIR}.ExampleQiimeMappingFile.txt
totalReads=`cat ../../Logs/${DIR}.phixOverall.stats.txt | grep "reads" | cut -f1 -d " ";`
phiXReads=`cat ../../Logs/${DIR}.phixOverall.stats.txt | grep "exactly" | cut -f5 -d " "`;
phiXReads1=`cat ../../Logs/${DIR}.phixOverall.stats.txt | grep ">1" | cut -f5 -d " "`;
phiXReads=$((${phiXReads}+${phiXReads1}));
barcoded=0;
for i in `cat ../../${DIR}.barcodeCounts.txt | cut -f2`; do barcoded=$(($barcoded + `echo $i`)); done;
strictMerged=`cat ../../Deliverables/${DIR}.Stats.Combined.txt | grep "Raw" -A2 | grep "Total" | cut -f2 -d " "`;
standardMerged=`bzcat ../../Deliverables/${DIR}.MergedRaw.fq.bz2 | grep "@HWI-" | wc -l`;
numSamples=`cat ../../Deliverables/${DIR}.Stats.Combined.txt | cut -f1 -d " " | wc -l`;
numSamples=$(($numSamples - 1));
totalReads=$((${totalReads}*250));
barcoded=$((${barcoded}*500));
phiXReads=$((${phiXReads}*250));
strictMerged=$((${strictMerged}*500));
standardMerged=$((${standardMerged}*500));
percentBarcoded=$((100*${barcoded}/${totalReads}));
percentPhiX=$((100*${phiXReads}/${totalReads}));
percentStrict=$((100*${strictMerged}/${barcoded}));
percentStandard=$((100*${standardMerged}/${barcoded}));
echo -e "PoolID\t#Samples\tMachine Reads\tBarcoded Reads\t%Barcoded Reads:Machine Reads\t#phiX Mapped\t%phiX:Machine Reads\tStrict Merged Reads\t%Strict Merged Reads:Barcoded Reads\tStandard Merged Reads\t%Standard Merged Reads:Barcoded Reads\n${DIR}\t${numSamples}\t${totalReads}\t${barcoded}\t${percentBarcoded}%\t${phiXReads}\t${percentPhiX}%\t${strictMerged}\t${percentStrict}%\t${standardMerged}\t${percentStandard}%" > ../../Deliverables/${DIR}.PoolStats.txt;
chmod -R 777 ../../Deliverables
