#!/bin/sh

THREADS=$1;

if [ -z "$THREADS" ];
	then THREADS=`grep -c ^processor /proc/cpuinfo`;
fi

mkdir ../../Deliverables;
DIR=`pwd -P | cut -f9 -d "/" | cut -f1 | sed 's/WorkDir//g'`;

cat ../SampleList | parallel -j$THREADS -I {} 'echo {}; usearch70 -fastq_mergepairs {}.1.fq -reverse {}.2.fq -fastq_minovlen 50 -fastq_maxdiffs 4 -fastq_truncqual 5 -fastqout {}.mergedRaw.fq';
cat ../SampleList | parallel -j$THREADS -I {} 'echo {}; usearch70 -fastq_filter {}.mergedRaw.fq -fastq_maxee .5 -fastqout {}.filteredRaw.fq';
cat ../SampleList | parallel -j$THREADS -I {} 'echo {}; usearch70 -fastq_mergepairs {}.1.fq -reverse {}.2.fq -fastq_minovlen 50 -fastq_maxdiffs 0 -fastq_minmergelen 252 -fastq_maxmergelen 254 -fastq_truncqual 5 -fastqout {}.merged.fq';
cat ../SampleList | parallel -j$THREADS -I {} 'echo {}; usearch70 -fastq_filter {}.merged.fq -fastq_maxee .05 -relabel "{}_" -fastqout {}.filtered.fq';
cat *.filteredRaw.fq > seqs.raw.fq;
cat *.filtered.fq > seqs.fq;
bowtie2 -x /users/mcwong/work/references/phix/phix -U seqs.raw.fq --end-to-end --very-sensitive --reorder -p $THREADS --un seqs.raw.filtered.fq -S /dev/null 2>../../Logs/phix.raw.prep.bleed.txt;
bowtie2 -x /users/mcwong/work/references/phix/phix -U seqs.fq --end-to-end --very-sensitive --reorder -p $THREADS --un seqs.filtered.fq -S /dev/null 2>../../Logs/phix.bleed.txt;
mv seqs.raw.filtered.fq MergedRaw.fq
mkdir ../split_libraries;
fq2fa seqs.filtered.fq seqs.fna;
mv seqs.fna ../split_libraries;
cat *.1.fq | sed 's/N:0:.*/N:0:/g' > Read1.fq;
cat *.2.fq | sed 's/N:0:.*/N:0:/g' > Read2.fq;
fq2fa MergedRaw.fq MergedRaw.fa;
cat Read1.fq Read2.fq > Reads.fq;
fq2fa Reads.fq Reads.fa;
usearch70 -fastq_mergepairs Read1.fq -reverse Read2.fq -fastaout temp.fa;
cat temp.fa >> Reads.fa
usearch70 -usearch_global MergedRaw.fa -db ~mcwong/jobs/Francisella/Francisella_V4.udb -strand plus -id .968 -uc FrancisellaMerged.uc -maxaccepts 0 -maxrejects 0 -threads $THREADS;
usearch70 -usearch_global Reads.fa -db ~mcwong/jobs/Francisella/Francisella_V4.udb -strand both -id .968 -uc FrancisellaRaw.uc -maxaccepts 0 -maxrejects 0 -threads $THREADS;
cat FrancisellaMerged.uc | cut -f9,10 | grep -v "*$" | cut -f1 | cut -f1 -d " " > Remove;
cat FrancisellaRaw.uc | cut -f9,10 | grep -v "*$" | cut -f1 | cut -f1 -d " " >> Remove;
cat MergedRaw.fq | perl /users/mcwong/fastqfilter.pl -v Remove > Temp.fq;
mv Temp.fq MergedRaw.fq;
cat Read1.fq | perl /users/mcwong/fastqfilter.pl -v Remove > Temp.fq;
mv Temp.fq Read1.fq;
cat Read2.fq | perl /users/mcwong/fastqfilter.pl -v Remove > Temp.fq;
mv Temp.fq Read2.fq;
pbzip2 -c Read1.fq > Read1.fq.bz2;
pbzip2 -c Read2.fq > Read2.fq.bz2;
pbzip2 -c MergedRaw.fq > MergedRaw.fq.bz2;
rm Read1.fq;
rm Read2.fq;
rm MergedRaw.fq;
cat ../SampleList | xargs -I {} mkdir {};
cat ../SampleList | parallel -I {} 'cat {}.1.fq | perl /users/mcwong/fastqfilter.pl -v Remove | pbzip2 -p1 -c > {}/{}.1.fq.bz2; cat {}.2.fq | perl /users/mcwong/fastqfilter.pl -v Remove | pbzip2 -p1 -c > {}/{}.2.fq.bz2';
rm seqs.raw.fq;
rm *.filteredRaw.fq;
rm *.merged.fq;
rm temp.fa
rm MergedRaw.fa;
rm Reads.fa;
rm Reads.fq;
rm FrancisellaMerged.uc;
rm FrancisellaRaw.uc;
rm seqs.filtered.fq;
rm seqs.fq;
rm *.filtered.fq;
~mcwong/gitRepo/16S_workflows/recoverBarcodesForRaw.pl Read1.fq.bz2 ../../${DIR}Barcodes/Project_${DIR}/Sample_${DIR}/${DIR}_NoIndex_L001_R2_001.fastq.gz | pbzip2 -c > RawReadsBarcodes.fq.bz2
~mcwong/gitRepo/16S_workflows/recoverBarcodesForRaw.pl MergedRaw.fq.bz2 ../../${DIR}Barcodes/Project_${DIR}/Sample_${DIR}/${DIR}_NoIndex_L001_R2_001.fastq.gz | pbzip2 -c > MergedBarcodes.fq.bz2
cd ..;
mkdir uparse;
usearch70 -derep_fulllength split_libraries/seqs.fna  -output uparse/derep.fna -sizeout -uc uparse/derep.uc;
usearch70 -sortbysize       uparse/derep.fna -output uparse/sorted.fa -minsize 2;
cp uparse/sorted.fa uparse/temp.fa
for i in {0.4,0.8,1.2,1.6,2.0,2.4,2.8,3.2};
do 
usearch70 -cluster_otus     uparse/temp.fa -otus   uparse/temp1.fa -otu_radius_pct $i -uc uparse/cluster_$i.uc -fastaout uparse/clustering.$i.fasta.out;
cat uparse/clustering.$i.fasta.out | grep "^>" | grep chimera | sed 's/^>//g' | sed -re 's/;n=.*up=/\t/g' | sed 's/;$//g' | tee -a uparse/chimeras.txt > uparse/chimeras.$i.txt;
cat uparse/clustering.$i.fasta.out | grep "^>" > uparse/uparseref.decisions.$i.txt;
rm uparse/clustering.$i.fasta.out;
mv uparse/temp1.fa uparse/temp.fa;
done
mv uparse/temp.fa uparse/otus1.fa
usearch70 -uchime_ref       uparse/otus1.fa  -db /users/mcwong/work/references/gold/gold.fa -strand plus -nonchimeras uparse/otus.fa -uchimeout uparse/uchimeref.uc;
usearch70 -usearch_global   uparse/otus.fa -db /stornext/snfs7/cmmr/dpsmith/silva/latest/silva_V4.udb -id .968 -strand plus -threads $THREADS -uc uparse/otus2taxa.uc -maxaccepts 0 -maxrejects 0;
cat uparse/derep.fna | grep -A1 "size=1;" | cut -f2 -d ">" | getSeq uparse/derep.fna > uparse/singletons.fna;
usearch70 -usearch_global uparse/singletons.fna -db uparse/sorted.fa -id .99 -uc uparse/singletons2otus.uc -strand plus -threads $THREADS -maxaccepts 32 -maxrejects 128 -minqt 1 -leftjust -rightjust -wordlength 12;
resolveIterativeUparse.pl uparse/cluster_*.uc uparse/singletons2otus.uc uparse/otus2taxa.uc --derep uparse/derep.uc --chimeras uparse/chimeras.txt --uchime uparse/uchimeref.uc --taxonomy /stornext/snfs7/cmmr/dpsmith/silva/latest/silva.map;
mv otu_table.biom uparse/;
biom summarize-table -i uparse/otu_table.biom -o uparse/stats.otu_table.txt;
cat split_libraries/seqs.fna | grep "^>" | cut -f1 -d "_" | cut -f2 -d ">" | sort | uniq -c > Stats.MergedReads.txt;
cat uparse/stats.otu_table.txt | tail -n +17 | sed 's/^ //g' | sed -re 's/: /\t/g' | sed 's/\.0$//g' > Stats.MappedReads.txt;
perl ~mcwong/gitRepo/16S_workflows/StatsComparisonMergedVsMapped.pl Stats.MergedReads.txt Stats.MappedReads.txt > Stats.Combined.txt
summarize_taxa.py -i uparse/otu_table.biom -o ../Deliverables/${DIR}.TaxaSummary
cp uparse/otu_table.biom ../Deliverables/${DIR}.otu_table.biom; 
cp Stats.Combined.txt ../Deliverables/${DIR}.Stats.Combined.txt;
cp SampleList ../Deliverables/${DIR}.SampleList;
cd Reads;
ls *.fq.bz2 | parallel -I {} mv {} ../../Deliverables/${DIR}.{}
cd ../..
cat samplesheet.${DIR}.csv | cut -f3,5 -d "," | tr "," "\t" > Deliverables/${DIR}.SampleSheet.txt
head -1 ~mcwong/IlluminaHeaderExample > Deliverables/${DIR}.ExampleQiimeMappingFile.txt
tail -n+2  Deliverables/${DIR}.SampleSheet.txt | sed -re 's/(.*)\t(.*)/\1\t\2\tGGACTACHVGGGTWTCTAAT\tGTGCCAGCMGCCGCGGTAA\t\1/g' >> Deliverables/${DIR}.ExampleQiimeMappingFile.txt
chmod -R 777 Deliverables
