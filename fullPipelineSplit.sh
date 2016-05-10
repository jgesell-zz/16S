#!/bin/sh

export READSDIR=$1;
export THREADS=$2;
CURRWORKDIR=`pwd`;

#If no thread count was passed, used all available threads on the cluster
if [ -z "${THREADS}" ];
	then export THREADS=`grep -c ^processor /proc/cpuinfo`;
fi

#If no reads directory was passed, sets it to the current directory
if [ -z "${READSDIR}" ];
	then export READSDIR=".";
fi

#Verify the existence of the temporary directory
if [ -d $(readlink -e $TMPDIR) ];
	then echo "Temporary Directory: ${TMPDIR}";
	else echo "Temporary Directory does not exist";
fi

#go into Reads dir
export READSDIR=`readlink -e $READSDIR`;
cd ${READSDIR};
export PROJECTID=`basename $(readlink -e ${READSDIR}/..) | sed 's/WorkDir//g'`;
echo ${PROJECTID};

## prepping the fastqs for processing ##
#decompress fastqs to temporary directory
cat  ${READSDIR}/../SampleList | parallel -j${THREADS} -I {} "bzcat {}.1.fq.bz2 | sed 's/1:N:0:.*/1:N:0:/g' > ${TMPDIR}/{}.1.fq";
cat  ${READSDIR}/../SampleList | parallel -j${THREADS} -I {} "bzcat {}.2.fq.bz2 | sed 's/2:N:0:.*/3:N:0:/g' > ${TMPDIR}/{}.2.fq";

#merge both standard merge and strict merge
cat  ${READSDIR}/../SampleList | parallel -j${THREADS} -I {} "echo '{} Standard Merge:'; usearch70 -fastq_mergepairs ${TMPDIR}/{}.1.fq -reverse ${TMPDIR}/{}.2.fq -fastq_minovlen 50 -fastq_maxdiffs 4 -fastq_truncqual 5 -fastqout ${TMPDIR}/{}.MergedStandard.fq; echo";
cat  ${READSDIR}/../SampleList | parallel -j${THREADS} -I {} "echo '{} Strict Merge:'; usearch70 -fastq_mergepairs ${TMPDIR}/{}.1.fq -reverse ${TMPDIR}/{}.2.fq -fastq_minovlen 50 -fastq_maxdiffs 0 -fastq_minmergelen 252 -fastq_maxmergelen 254 -fastq_truncqual 5 -fastqout ${TMPDIR}/{}.MergedStrict.fq; echo";

#filter three  ways (raw merged, standard merged, strict merged)
cat  ${READSDIR}/../SampleList | parallel -j${THREADS} -I {} "echo '{} Filter Raw Merge:'; usearch70 -fastq_filter ${TMPDIR}/{}.MergedStandard.fq -fastq_maxee .5 -fastqout ${TMPDIR}/{}.filteredRaw.fq; echo";
cat  ${READSDIR}/../SampleList | parallel -j${THREADS} -I {} "echo '{} Filtered Standard Merge:'; usearch70 -fastq_filter ${TMPDIR}/{}.MergedStandard.fq -fastq_maxee .5 -relabel \"{}_\" -fastqout ${TMPDIR}/{}.filteredStandard.fq; echo";
cat  ${READSDIR}/../SampleList | parallel -j${THREADS} -I {} "echo '{} Filter Strict Merge:'; usearch70 -fastq_filter ${TMPDIR}/{}.MergedStrict.fq -fastq_maxee .05 -relabel \"{}_\" -fastqout ${TMPDIR}/{}.filteredStrict.fq; echo";

#concatenate demultiplexed fastqs into monolithic variants
cat ${TMPDIR}/*.filteredStrict.fq > ${TMPDIR}/seqs.strict.fq &
cat ${TMPDIR}/*.filteredRaw.fq > ${TMPDIR}/seqs.raw.fq &
cat ${TMPDIR}/*.filteredStandard.fq > ${TMPDIR}/seqs.standard.fq &
wait;

#filter out phiX bleed
bowtie2 -x ${PHIXDB} -U ${TMPDIR}/seqs.strict.fq --end-to-end --very-sensitive --reorder -p ${THREADS} --un ${TMPDIR}/seqs.strict.filtered.fq -S /dev/null 2>&1;
bowtie2 -x ${PHIXDB} -U ${TMPDIR}/seqs.raw.fq --end-to-end --very-sensitive --reorder -p ${THREADS} --un ${TMPDIR}/seqs.raw.filtered.fq -S /dev/null 2>&1;
bowtie2 -x ${PHIXDB} -U ${TMPDIR}/seqs.standard.fq --end-to-end --very-sensitive --reorder -p ${THREADS} --un ${TMPDIR}/seqs.standard.filtered.fq -S /dev/null 2>&1;


#construct the fastas for uparse
mkdir  ${READSDIR}/../split_libraries;
fq2fa ${TMPDIR}/seqs.strict.filtered.fq  ${READSDIR}/../split_libraries/Strict.seqs.fna &
fq2fa ${TMPDIR}/seqs.standard.filtered.fq  ${READSDIR}/../split_libraries/Standard.seqs.fna &
wait;

#simultaneous uparse
for j in {Strict,Standard}; do echo $j; done | parallel -I {} '
mkdir ${TMPDIR}/uparse{};
usearch70 -derep_fulllength  ${READSDIR}/../split_libraries/{}.seqs.fna  -output ${TMPDIR}/uparse{}/derep.fna -sizeout -uc ${TMPDIR}/uparse{}/derep.uc 2>&1;
usearch70 -sortbysize ${TMPDIR}/uparse{}/derep.fna -output ${TMPDIR}/uparse{}/sorted.fa -minsize 2;
cp ${TMPDIR}/uparse{}/sorted.fa ${TMPDIR}/uparse{}/temp.fa;
for i in {0.4,0.8,1.2,1.6,2.0,2.4,2.8,3.2};
	do usearch70 -cluster_otus ${TMPDIR}/uparse{}/temp.fa -otus   ${TMPDIR}/uparse{}/temp1.fa -otu_radius_pct $i -uc ${TMPDIR}/uparse{}/cluster_$i.uc -fastaout ${TMPDIR}/uparse{}/clustering.$i.fasta.out;
	cat ${TMPDIR}/uparse{}/clustering.$i.fasta.out | grep "^>" | grep chimera | sed "s/^>//g" | sed -re "s/;n=.*up=/\t/g" | sed "s/;$//g" | tee -a ${TMPDIR}/uparse{}/chimeras.txt > ${TMPDIR}/uparse{}/chimeras.$i.txt;
	cat ${TMPDIR}/uparse{}/clustering.$i.fasta.out | grep "^>" > ${TMPDIR}/uparse{}/uparse{}ref.decisions.$i.txt;
	rm ${TMPDIR}/uparse{}/clustering.$i.fasta.out;
	mv ${TMPDIR}/uparse{}/temp1.fa ${TMPDIR}/uparse{}/temp.fa;
done;
mv ${TMPDIR}/uparse{}/temp.fa ${TMPDIR}/uparse{}/otus1.fa;
usearch70 -uchime_ref ${TMPDIR}/uparse{}/otus1.fa  -db ${GOLD} -strand plus -nonchimeras ${TMPDIR}/uparse{}/otus.fa -uchimeout ${TMPDIR}/uparse{}/uchimeref.uc;
usearch70 -usearch_global ${TMPDIR}/uparse{}/otus.fa -db ${SILVA}/silva_V4.udb -id .968 -strand plus -threads ${THREADS} -uc ${TMPDIR}/uparse{}/otus2taxa.uc -maxaccepts 0 -maxrejects 0;
cat ${TMPDIR}/uparse{}/derep.fna | grep -A1 "size=1;" | cut -f2 -d ">" | ${GITREPO}/Miscellaneous/getSeq ${TMPDIR}/uparse{}/derep.fna > ${TMPDIR}/uparse{}/singletons.fna;
usearch70 -usearch_global ${TMPDIR}/uparse{}/singletons.fna -db ${TMPDIR}/uparse{}/sorted.fa -id .99 -uc ${TMPDIR}/uparse{}/singletons2otus.uc -strand plus -threads ${THREADS} -maxaccepts 32 -maxrejects 128 -minqt 1 -leftjust -rightjust -wordlength 12;
cd ${TMPDIR}/uparse{};
${WONGGITREPO}/16S_workflows/resolveIterativeUparse.pl ${TMPDIR}/uparse{}/cluster_*.uc ${TMPDIR}/uparse{}/singletons2otus.uc ${TMPDIR}/uparse{}/otus2taxa.uc --derep ${TMPDIR}/uparse{}/derep.uc --chimeras ${TMPDIR}/uparse{}/chimeras.txt --uchime ${TMPDIR}/uparse{}/uchimeref.uc --taxonomy ${SILVA}/silva.map;
biom summarize-table -i ${TMPDIR}/uparse{}/otu_table.biom -o ${TMPDIR}/uparse{}/Stats.{}Merge.otu_table.txt;
if [ $? -eq 0 ];
	then echo "biom {} Stats run successful";
	else echo "biom {} Stats run failed";
	exit 1;
fi;
var=$(expr `cat ${TMPDIR}/uparse{}/Stats.{}Merge.otu_table.txt | grep -n "Counts/sample detail" | cut -f1 -d ":"` + 1);
if [ $var -eq 1 ];
	then echo "No start-line found, terminating";
	exit 1;
fi;
cat ${TMPDIR}/uparse{}/Stats.{}Merge.otu_table.txt | tail -n +$var | sed "s/^ //g" | sed -re "s/: /\t/g" | sed "s/\.0$//g" >  ${READSDIR}/../Stats.{}Merge.MappedReads.txt;
cat  ${READSDIR}/../split_libraries/{}.seqs.fna | grep "^>" | cut -f1 -d "_" | cut -f2 -d ">" | sort | uniq -c >  ${READSDIR}/../Stats.{}Merge.MergedReads.txt;
perl ${WONGGITREPO}/16S_workflows/StatsComparisonMergedVsMapped.pl  ${READSDIR}/../Stats.{}Merge.MergedReads.txt  ${READSDIR}/../Stats.{}Merge.MappedReads.txt >  ${READSDIR}/../Stats.{}Merge.Combined.txt;
tar -cvjf  ${READSDIR}/../uparse{}.tar.bz2 ${TMPDIR}/uparse{};
' &
bigJob=`jobs -p`;

## construct the deliverables ##
#set up francisella filter
mv ${TMPDIR}/seqs.raw.filtered.fq ${TMPDIR}/MergedStandard.fq;
cat  ${READSDIR}/../SampleList | parallel -j1 -I {} "cat ${TMPDIR}/{}.1.fq" > ${TMPDIR}/Read1.fq;
cat  ${READSDIR}/../SampleList | parallel -j1 -I {} "cat ${TMPDIR}/{}.2.fq" > ${TMPDIR}/Read2.fq;
fq2fa ${TMPDIR}/MergedStandard.fq ${TMPDIR}/MergedStandard.fa;
cat ${TMPDIR}/Read1.fq ${TMPDIR}/Read2.fq > ${TMPDIR}/Reads.fq;
fq2fa ${TMPDIR}/Reads.fq ${TMPDIR}/Reads.fa;
usearch70 -fastq_mergepairs ${TMPDIR}/Read1.fq -reverse ${TMPDIR}/Read2.fq -fastaout ${TMPDIR}/temp.fa;
cat ${TMPDIR}/temp.fa >> ${TMPDIR}/Reads.fa;
cat ${TMPDIR}/MergedStandard.fa >> ${TMPDIR}/Reads.fa;

#filter francisella
usearch70 -usearch_global ${TMPDIR}/Reads.fa -db ${FRANCISELLA}/Francisella_V4.udb -strand both -id .968 -uc ${TMPDIR}/Francisella.uc -maxaccepts 0 -maxrejects 0 -threads ${THREADS};

#remove francisella
cat ${TMPDIR}/Francisella.uc | cut -f9,10 | grep -v "*$" | cut -f1 | cut -f1 -d " " > ${TMPDIR}/Remove;
cat ${TMPDIR}/MergedStandard.fq | perl ${GITREPO}/Miscellaneous/fastqfilter.pl -v ${TMPDIR}/Remove > ${TMPDIR}/TempMerged.fq &
cat ${TMPDIR}/Read1.fq | perl ${GITREPO}/Miscellaneous/fastqfilter.pl -v ${TMPDIR}/Remove > ${TMPDIR}/Temp1.fq &
cat ${TMPDIR}/Read2.fq | perl ${GITREPO}/Miscellaneous/fastqfilter.pl -v ${TMPDIR}/Remove > ${TMPDIR}/Temp2.fq &
for i in `jobs -p | grep -v $bigJob`; do pidlist=`echo -e "\`echo $pidlist\`\n\`echo  $((i + 1))\`"`; done
wait $pidlist

#Move the temporary files to their final versions
mv ${TMPDIR}/TempMerged.fq ${TMPDIR}/MergedStandard.fq;
mv ${TMPDIR}/Temp1.fq ${TMPDIR}/Read1.fq;
mv ${TMPDIR}/Temp2.fq ${TMPDIR}/Read2.fq;

#compress  deliverables
pbzip2 -p${THREADS} ${TMPDIR}/Read1.fq
pbzip2 -p${THREADS} ${TMPDIR}/Read2.fq
pbzip2 -p${THREADS} ${TMPDIR}/MergedStandard.fq
cat  ${READSDIR}/../SampleList | xargs -I {} mkdir {};
cat  ${READSDIR}/../SampleList | parallel -j${THREADS} -I {} 'cat ${TMPDIR}/{}.2.fq | perl ${GITREPO}/Miscellaneous/fastqfilter.pl -v ${TMPDIR}/Remove | pbzip2 -p1 -c > {}/{}.2.fq.bz2';
cat  ${READSDIR}/../SampleList | parallel -j${THREADS} -I {} 'cat ${TMPDIR}/{}.1.fq | perl ${GITREPO}/Miscellaneous/fastqfilter.pl -v ${TMPDIR}/Remove | pbzip2 -p1 -c > {}/{}.1.fq.bz2';

#recover barcodes for deliverables
${WONGGITREPO}/16S_workflows/recoverBarcodesForRaw.pl ${TMPDIR}/Read1.fq.bz2  ${READSDIR}/../../${PROJECTID}Barcodes/Project_${PROJECTID}/Sample_${PROJECTID}/${PROJECTID}_NoIndex_L001_R2_001.fastq.gz | pbzip2 -p${THREADS} -c > ${TMPDIR}/RawReadsBarcodes.fq.bz2;
${WONGGITREPO}/16S_workflows/recoverBarcodesForRaw.pl ${TMPDIR}/MergedStandard.fq.bz2   ${READSDIR}/../../${PROJECTID}Barcodes/Project_${PROJECTID}/Sample_${PROJECTID}/${PROJECTID}_NoIndex_L001_R2_001.fastq.gz | pbzip2 -p${THREADS} -c > ${TMPDIR}/MergedBarcodes.fq.bz2;

mkdir  ${READSDIR}/../../Deliverables;

#wait for otu_table construction to finish
wait; 

#move files to their destination for further analysis
for i in `ls ${TMPDIR}/*.fq.bz2`; do name=`basename $i`; cp $i  ${READSDIR}/../../Deliverables/${PROJECTID}.${name};done;
cp ${TMPDIR}/uparseStandard/otu_table.biom  ${READSDIR}/../../Deliverables/${PROJECTID}.Standard.otu_table.biom; 
cp ${TMPDIR}/uparseStrict/otu_table.biom  ${READSDIR}/../../Deliverables/${PROJECTID}.Strict.otu_table.biom; 
cp  ${READSDIR}/../Stats.StrictMerge.Combined.txt  ${READSDIR}/../../Deliverables/${PROJECTID}.Stats.StrictMerge.Combined.txt;
cp  ${READSDIR}/../Stats.StandardMerge.Combined.txt  ${READSDIR}/../../Deliverables/${PROJECTID}.Stats.StandardMerge.Combined.txt;
cp  ${READSDIR}/../SampleList  ${READSDIR}/../../Deliverables/${PROJECTID}.SampleList.txt;
cat  ${READSDIR}/../../samplesheet.${PROJECTID}.csv | grep -f ${READSDIR}/../SampleList | cut -f3,5 -d "," | tr "," "\t" >  ${READSDIR}/../../Deliverables/${PROJECTID}.SampleSheet.txt;
head -1 ${GITREPO}/Miscellaneous/IlluminaHeaderExample >  ${READSDIR}/../../Deliverables/${PROJECTID}.ExampleQiimeMappingFile.txt;
tail -n+1  ${READSDIR}/../../Deliverables/${PROJECTID}.SampleSheet.txt | sed -re 's/(.*)\t(.*)/\1\t\2\tGGACTACHVGGGTWTCTAAT\tGTGCCAGCMGCCGCGGTAA\t\1/g' >>  ${READSDIR}/../../Deliverables/${PROJECTID}.ExampleQiimeMappingFile.txt;
chmod -R 777  ${READSDIR}/../../Deliverables;

#return to working directory when script was launched
cd $CURRWORKDIR;

#exit without error status once completed
exit 0;
