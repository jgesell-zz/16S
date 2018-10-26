#!/bin/sh

umask 002;
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
if [ -d $(readlink -e ${TMPDIR}) ];
then echo "Temporary Directory: ${TMPDIR}";
else echo "Temporary Directory does not exist";
fi

#go into Reads dir
export READSDIR=`readlink -e ${READSDIR}`;
cd ${READSDIR};
export PROJECTID=`basename $(readlink -e ${READSDIR}/..) | sed 's/WorkDir//g'`;
echo "Current Reads Directory: ${READSDIR}";
echo "Current Project ID: ${PROJECTID}";
THREADSPLIT=$[THREADS / 2];

## prepping the fastqs for processing ##
#decompress fastqs to temporary directory
cat ${READSDIR}/../SampleList | parallel -j${THREADS} -I {} "bzcat {}.1.fq.bz2 | sed 's/1:N:0:.*/1:N:0:/g' > ${TMPDIR}/{}.1.fq";
cat ${READSDIR}/../SampleList | parallel -j${THREADS} -I {} "bzcat {}.2.fq.bz2 | sed 's/2:N:0:.*/3:N:0:/g' > ${TMPDIR}/{}.2.fq";

#Verify that the barcode counts file exists, and if it does not, create it
if [ ! -z ${READSDIR}/../../${PROJECTID}.barcodeCounts.txt ];
then for i in `cat ${READSDIR}/../SampleList`;
do count=`cat ${TMPDIR}/${i}.1.fq | wc -l`;
count=$[count / 4];
echo -e "${i}\t${count}";
done > ${READSDIR}/../../${PROJECTID}.barcodeCounts.txt;
fi; 

#filter three ways (raw merged, standard merged, strict merged)
cat ${READSDIR}/../SampleList | parallel -j${THREADS} -I {} "echo '{} Filter Raw Merge'; usearch81 -fastq_filter ${TMPDIR}/{}.1.fq -fastq_maxee_rate .05 -threads 1 -fastqout ${TMPDIR}/{}.filteredRaw.fq && rm ${TMPDIR}/{}.MergedStandard.fq; echo";
cat ${READSDIR}/../SampleList | parallel -j${THREADS} -I {} "echo '{} Filter Strict Merge'; usearch81 -fastq_filter ${TMPDIR}/{}.1.fq -fastq_maxee_rate .05 -threads 1 -relabel \"{}_\" -fastqout ${TMPDIR}/{}.filteredStrict.fq && rm ${TMPDIR}/{}.Merged_Reads.fq; echo";

#concatenate demultiplexed fastqs into monolithic variants
cat ${TMPDIR}/*.filteredStrict.fq > ${TMPDIR}/seqs.strict.fq &
cat ${TMPDIR}/*.filteredRaw.fq > ${TMPDIR}/seqs.raw.fq &
cat ${READSDIR}/../../${PROJECTID}.barcodeCounts.txt | grep -f ${READSDIR}/../SampleList > ${TMPDIR}/${PROJECTID}.barcodeCounts.txt &
wait;
rm ${TMPDIR}/*.filteredStrict.fq ${TMPDIR}/*.filteredRaw.fq &

#filter out phiX bleed
bowtie2 -x ${PHIXDB} -U ${TMPDIR}/seqs.strict.fq --end-to-end --very-sensitive --reorder -p ${THREADSPLIT} --un ${TMPDIR}/seqs.strict.filtered.fq -S /dev/null 2>&1 && rm ${TMPDIR}/seqs.strict.fq &
bowtie2 -x ${PHIXDB} -U ${TMPDIR}/seqs.raw.fq --end-to-end --very-sensitive --reorder -p ${THREADSPLIT} --un ${TMPDIR}/seqs.raw.filtered.fq -S /dev/null 2>&1 && rm ${TMPDIR}/seqs.raw.fq &
wait;

#construct the fastas for uparse
mkdir ${READSDIR}/../split_libraries;
mkdir -p ${READSDIR}/../../Deliverables.Read1Only;
fq2fa ${TMPDIR}/seqs.strict.filtered.fq ${READSDIR}/../split_libraries/Strict.seqs.fna;

#simultaneous uparse
echo "Strict" | parallel -I {} '
mkdir ${TMPDIR}/uparse{};
usearch70 -derep_fulllength ${READSDIR}/../split_libraries/{}.seqs.fna -output ${TMPDIR}/uparse{}/derep.fna -sizeout -uc ${TMPDIR}/uparse{}/derep.uc 2>&1;
usearch70 -sortbysize ${TMPDIR}/uparse{}/derep.fna -output ${TMPDIR}/uparse{}/sorted.fa -minsize 2;
cp ${TMPDIR}/uparse{}/sorted.fa ${TMPDIR}/uparse{}/temp.fa;
for i in {0.4,0.8,1.2,1.6,2.0,2.4,2.8,3.2};
do usearch70 -cluster_otus ${TMPDIR}/uparse{}/temp.fa -otus ${TMPDIR}/uparse{}/temp1.fa -otu_radius_pct $i -uc ${TMPDIR}/uparse{}/cluster_$i.uc -fastaout ${TMPDIR}/uparse{}/clustering.$i.fasta.out;
cat ${TMPDIR}/uparse{}/clustering.$i.fasta.out | grep "^>" | grep chimera | sed "s/^>//g" | sed -re "s/;n=.*up=/\t/g" | sed "s/;$//g" | tee -a ${TMPDIR}/uparse{}/chimeras.txt > ${TMPDIR}/uparse{}/chimeras.$i.txt;
cat ${TMPDIR}/uparse{}/clustering.$i.fasta.out | grep "^>" > ${TMPDIR}/uparse{}/uparse{}ref.decisions.$i.txt;
rm ${TMPDIR}/uparse{}/clustering.$i.fasta.out;
mv ${TMPDIR}/uparse{}/temp1.fa ${TMPDIR}/uparse{}/temp.fa;
done;
mv ${TMPDIR}/uparse{}/temp.fa ${TMPDIR}/uparse{}/otus1.fa;
usearch70 -uchime_ref ${TMPDIR}/uparse{}/otus1.fa -db ${GOLD} -strand plus -uchimeout ${TMPDIR}/uparse{}/uchimeref.uc;
usearch70 -usearch_global ${TMPDIR}/uparse{}/otus.fa -db ${SILVA}/silva.udb -id .968 -strand plus -threads ${THREADS} -uc ${TMPDIR}/uparse{}/otus2taxa.uc -maxaccepts 0 -maxrejects 0;
cat ${TMPDIR}/uparse{}/derep.fna | grep -A1 "size=1;" | cut -f2 -d ">" | ${GITREPO}/Miscellaneous/getSeq ${TMPDIR}/uparse{}/derep.fna > ${TMPDIR}/uparse{}/singletons.fna;
usearch70 -usearch_global ${TMPDIR}/uparse{}/singletons.fna -db ${TMPDIR}/uparse{}/sorted.fa -id .99 -uc ${TMPDIR}/uparse{}/singletons2otus.uc -strand plus -threads ${THREADS} -maxaccepts 32 -maxrejects 128 -minqt 1 -leftjust -rightjust -wordlength 12;
cd ${TMPDIR}/uparse{};
${GITREPO}/Miscellaneous/resolveIterativeUparse.pl ${TMPDIR}/uparse{}/cluster_*.uc ${TMPDIR}/uparse{}/singletons2otus.uc ${TMPDIR}/uparse{}/otus2taxa.uc --derep ${TMPDIR}/uparse{}/derep.uc --chimeras ${TMPDIR}/uparse{}/chimeras.txt --uchime ${TMPDIR}/uparse{}/uchimeref.uc --taxonomy ${SILVA}/silva.map;
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
cat ${TMPDIR}/uparse{}/Stats.{}Merge.otu_table.txt | tail -n +$var | sed "s/^ //g" | sed -re "s/: /\t/g" | sed "s/\.0$//g" > ${READSDIR}/../Stats.{}Merge.MappedReads.txt &
cat ${READSDIR}/../split_libraries/{}.seqs.fna | grep "^>" | cut -f1 -d "_" | cut -f2 -d ">" | sort | uniq -c > ${READSDIR}/../Stats.{}Merge.MergedReads.txt &
wait;
perl ${GITREPO}/16S/StatsComparisonMergedVsMapped.pl ${TMPDIR}/${PROJECTID}.barcodeCounts.txt ${READSDIR}/../Stats.{}Merge.MergedReads.txt ${READSDIR}/../Stats.{}Merge.MappedReads.txt > ${READSDIR}/../Stats.{}Merge.Combined.txt;
for i in `cat ${TMPDIR}/uparse{}/otus.fa | perl -pe "s/;\n/;/g"`; do sample=`echo ${i} | cut -f1 -d ";" | sed -e "s:^>::g"`; count=`echo ${i} | cut -f2 -d ";"`; otu=`cat ${TMPDIR}/uparse{}/reads2otus.txt | grep -w ${sample} | cut -f2`; seq=`echo ${i} | cut -f3 -d ";"`; if  grep -q -w "${otu}" otu_table.biom; then echo -e ">${sample};${otu};${count};\n${seq}"; fi; done | pbzip2 -c > ${READSDIR}/../../Deliverables.Read1Only/CentroidInformation.fa.bz2 &
cp ${TMPDIR}/uparse{}/otu_table.biom ${READSDIR}/../../Deliverables.Read1Only/OTU_Table.biom &
cp ${READSDIR}/../Stats.{}Merge.Combined.txt ${READSDIR}/../../Deliverables.Read1Only/Read_QC.txt &
wait;
tar -cvf ${READSDIR}/../uparse{}.tar.bz2 -C ${TMPDIR} uparse{} -I pbzip2 && rm -rf ${TMPDIR}/uparse{};
' &
bigJob=`jobs -p`;

## construct the deliverables ##
#set up francisella filter
mv ${TMPDIR}/seqs.strict.filtered.fq ${TMPDIR}/Merged_Reads.fq;
mv ${TMPDIR}/seqs.raw.filtered.fq ${TMPDIR}/Raw_Merged_Reads.fq;
fq2fa ${TMPDIR}/Raw_Merged_Reads.fq ${TMPDIR}/Raw_Merged_Reads.fa;
cat ${READSDIR}/../SampleList | parallel -j1 -I {} "cat ${TMPDIR}/{}.1.fq" > ${TMPDIR}/Raw_Read1.fq;
cat ${READSDIR}/../SampleList | parallel -j1 -I {} "cat ${TMPDIR}/{}.2.fq" > ${TMPDIR}/Raw_Read3.fq;
cat ${TMPDIR}/Raw_Read1.fq ${TMPDIR}/Raw_Read3.fq > ${TMPDIR}/Raw_Reads.fq;
fq2fa ${TMPDIR}/Raw_Reads.fq ${TMPDIR}/Raw_Reads.fa;
usearch70 -fastq_mergepairs ${TMPDIR}/Raw_Read1.fq -reverse ${TMPDIR}/Raw_Read3.fq -fastaout ${TMPDIR}/temp.fa;
cat ${TMPDIR}/temp.fa >> ${TMPDIR}/Raw_Reads.fa;
cat ${TMPDIR}/Raw_Merged_Reads.fa >> ${TMPDIR}/Raw_Reads.fa;

#Move the raw reads into Deliverables
cat ${READSDIR}/../SampleList | xargs -I {} mkdir -p ${READSDIR}/../../Deliverables.Read1Only/RawSequences/{};
cat ${READSDIR}/../SampleList | parallel -j${THREADS} -I {} 'cat ${TMPDIR}/{}.1.fq | pbzip2 -p1 -c > ${READSDIR}/../../Deliverables.Read1Only/RawSequences/{}/{}.1.fq.bz2 && rm ${TMPDIR}/{}.1.fq';
cat ${READSDIR}/../SampleList | parallel -j${THREADS} -I {} 'cat ${TMPDIR}/{}.2.fq | pbzip2 -p1 -c > ${READSDIR}/../../Deliverables.Read1Only/RawSequences/{}/{}.2.fq.bz2 && rm ${TMPDIR}/{}.2.fq';


#Move the temporary files to their final versions
#compress deliverables
pbzip2 -f -p${THREADS} ${TMPDIR}/Raw_Read1.fq;
pbzip2 -f -p${THREADS} ${TMPDIR}/Raw_Read3.fq;
pbzip2 -f -p${THREADS} ${TMPDIR}/Merged_Reads.fq;

#recover barcodes for deliverables
${GITREPO}/Miscellaneous/recoverBarcodesForRaw.pl ${TMPDIR}/Raw_Read1.fq.bz2 ${READSDIR}/../../${PROJECTID}Barcodes/Project_${PROJECTID}/Sample_${PROJECTID}/${PROJECTID}_NoIndex_L001_R2_001.fastq.gz | pbzip2 -p${THREADSPLIT} -c > ${TMPDIR}/Raw_Read2_Barcodes.fq.bz2 &
${GITREPO}/Miscellaneous/recoverBarcodesForRaw.pl ${TMPDIR}/Merged_Reads.fq.bz2 ${READSDIR}/../../${PROJECTID}Barcodes/Project_${PROJECTID}/Sample_${PROJECTID}/${PROJECTID}_NoIndex_L001_R2_001.fastq.gz | pbzip2 -p${THREADSPLIT} -c > ${TMPDIR}/Merged_Barcodes.fq.bz2 &
wait;
#move files to their destination for further analysis
cp ${TMPDIR}/*.fq.bz2 ${READSDIR}/../../Deliverables.Read1Only/ && rm ${TMPDIR}/*.fq.bz2;
head -1 ${GITREPO}/Miscellaneous/IlluminaHeaderExample > ${READSDIR}/../../Deliverables.Read1Only/Demultiplex_Sheet.txt;
cat  ${READSDIR}/../../samplesheet.${PROJECTID}.csv | grep -f ${READSDIR}/../SampleList | cut -f3,5 -d "," | tr "," "\t" | tail -n+1 | sed -re 's/(.*)\t(.*)/\1\t\2\tGGACTACHVGGGTWTCTAAT\tGTGCCAGCMGCCGCGGTAA\t\1/g' >> ${READSDIR}/../../Deliverables.Read1Only/Demultiplex_Sheet.txt;
cp ${GITREPO}/16S/CMMR16SV4Pipeline.md ${READSDIR}/../../Deliverables.Read1Only/;
biom convert -i ${READSDIR}/../../Deliverables.Read1Only/OTU_Table.biom -o ${READSDIR}/../../Deliverables.Read1Only/OTU_Table.tsv --to-tsv --header-key taxonomy;
/cmmr/bin/Rscript /cmmr/bin/deliver_folder.r -f ${READSDIR}/../../Deliverables.Read1Only -t ${SILVA}/silva.tre -n ${THREADS};
chmod -R 755 ${READSDIR}/../../Deliverables.Read1Only;
if [ -r "${READSDIR}/../../Deliverables.Read1Only/ProjectData.rds" ];
then collab=`readlink -e ${READSDIR} | cut -f5 -d "/"`;
pool=`readlink -e ${READSDIR} | cut -f6 -d "/"`;
if [ "${collab}" != "StatsProject" ];
then echo -e "${collab} ${pool} has completed running thru the 16S V4 pipeline.  Attached are the read statistics for this run.\nAll other deliverables can be found on the CMMR cluster at the following location:\t`readlink -e ${READSDIR}/../../Deliverables.Read1Only`" | mail -a ${READSDIR}/../../Deliverables.Read1Only/Read_QC.txt -s "${collab} ${pool} has completed" gesell@bcm.edu,mcross@bcm.edu,Nadim.Ajami@bcm.edu,carmical@bcm.edu;
elif [ "${collab}" = "StatsProject" ];
then echo -e "${collab} ${pool} has completed running thru the 16S V4 pipeline.  Attached are the read statistics for this run.\nAll other deliverables can be found on the CMMR cluster at the following location:\t`readlink -e ${READSDIR}/../../Deliverables.Read1Only`" | mail -a ${READSDIR}/../../Deliverables.Read1Only/Read_QC.txt -s "${collab} ${pool} has completed" gesell@bcm.edu, carmical@bcm.edu;
fi;
else echo -e "${collab} ${pool} run failed, please check reason" | mail -a ${READSDIR}/../../Deliverables.Read1Only/Read_QC.txt -s "${collab} ${pool} has failed" gesell@bcm.edu;
fi;

#return to working directory when script was launched
cd ${CURRWORKDIR};

#exit without error status once completed
exit 0;
