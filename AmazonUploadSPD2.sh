#!/bin/sh
export time=$1;
export link=$2

if [ -z "${time}" ];
then export time="30";
fi;

if [ -z "${link}" ];
then export link=`pwd -P`;
fi;

if [ -z "${TMPDIR}" ];
then echo "Error: Please run this program thru the queueing system!";
exit 1;
fi;

cd ${link};
file=`pwd -P | cut -f5,6 -d "/" | tr "/" "."`.SPD2
md5sum CMMR16SV4Pipeline.md Demultiplex_Sheet.txt Merged_Barcodes.fq.bz2 Merged_Reads.fq.bz2 Raw_Read1.fq.bz2 Raw_Read2_Barcodes.fq.bz2 Raw_Read3.fq.bz2 OTU_Table.biom Read_QC.txt CentroidInformation.fa.bz2 RawSequences/* Merged_Barcodes.fq.bz2 Merged_Reads.fq.bz2 OTU_Table.tsv OTU_Table.tre Summary_Tables.xlsx > SPD2.Upload.MD5Sums.txt;
zip -0r ${TMPDIR}/${file}.zip CMMR16SV4Pipeline.md Demultiplex_Sheet.txt Merged_Barcodes.fq.bz2 Merged_Reads.fq.bz2 Raw_Read1.fq.bz2 Raw_Read2_Barcodes.fq.bz2 Raw_Read3.fq.bz2 OTU_Table.biom Read_QC.txt CentroidInformation.fa.bz2 RawSequences/* Merged_Barcodes.fq.bz2 Merged_Reads.fq.bz2 OTU_Table.tsv OTU_Table.tre Summary_Tables.xlsx SPD2.Upload.MD5Sums.txt;
sum=`md5sum ${TMPDIR}/${file}.zip`;
aws s3 cp ${TMPDIR}/${file}.zip s3://jplab/share/${time}d/;
if [ $? -ne 0 ];
then echo -e "${file} failed to upload to S3, please check error logs for the reason why." | mail -s "${file} Upload Failure" ${USER}@bcm.edu;
else link=`aws s3 presign s3://jplab/share/${time}d/${file}.zip --expires-in $[time * 24 * 60 * 60]`;
echo -e "File Name\tMD5Sum\tLink\n${file}\t${sum}\t${link}" > SPD2.Upload.Link.txt;
echo -e "${file} successfully uploaded to the Amazon Cloud.\nLink:\t${link}\nMD5Sum:\t${sum}" | mail -s "${file} Upload Complete" ${USER}@bcm.edu;
fi;
exit 0;
