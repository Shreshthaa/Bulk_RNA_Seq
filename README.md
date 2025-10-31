# Bulk_RNA_Seq
This Project has been produced to demonstrate A simple yet reproducible bulk RNA seq pipeline.
Project title - Bulk_RMA_seq_Analysis
Created by - Shreshtha Shukla.  
 Database - The dataset used was taken from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158550) (accession id-GSE158550).Reference paper from [Nature communications](https://www.nature.com/articles/s41467-018-08133-6#content).

 ## Project Workflow
*This repository documents a beginner-friendly, reproducible Bulk RNA-Seq analysis workflow (download → QC → trimming → alignment → QC → counting → counts matrix preparation). The repository stores commands exactly as executed by the author, brief explanations, small result files (counts, QC HTMLs), scripts (text) and screenshots of terminal outputs so learners can reproduce the pipeline.*
## Step by step process 
Our goal is to find differentially expressed genes in response to hypoxia for the LNCaP cell lines. original database includes LNCap and PC3 cell lines for normoxia and hupoxia conditions but due to minimal storage, we have performed it on LNCap cell line exclusively. The entire pipeline has been divided into two major parts which is preprocessing and QC of data and expression analysis using R studio.
## Preprocessing, extraction of data and Quality control raw files.
# Installing SRA toolkit and fetching the SRR files.
```bash
sudo apt install sra-toolkit
```
```bash
prefetch SRR7179504
```
# Converting the .sra file into .fastq and further saving the output into fastq folder. The files are gzipped for efficienty using the storage.
```bash
fastq-dump --outdir fastq --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip ~/SRR7179504/SRR7179504.sra
```
# Creating an automated python script for converting all the .sra files into .fastq gzips
the specific code does more than that, It concatenates the SRR files into combined samples and renames the .fastq gzipped files into more meaningfull abbreviations.
```bash
python3 fastq_download.py
```
<img width="1054" height="46" alt="image" src="https://github.com/user-attachments/assets/cd1ad1fd-d842-4b8b-a1de-2722e5515c99" />

# Fastqc and Multiqc
in this step we move further to do a quality control assessment of the sequencing data we are working on : GC content,  understanding the whisker plots and number of errors per read etc. 
```bash
sudo apt install fastqc
```
```bash
fastqc fastq/*.fastq.gz -o fastqc_results/ --threads 8
```
The results are as follows:
<img width="798" height="485" alt="Screenshot 2025-10-31 161103" src="https://github.com/user-attachments/assets/8f9e538e-23c9-4834-9bad-acd6859043d3" />
a basic information chat indication towards the sequencer information( illumina hiseq 200) and read lenght 76 indicating a single end sequencing.
<img width="1315" height="675" alt="Screenshot 2025-10-31 161131" src="https://github.com/user-attachments/assets/3c7803d9-2d60-44b4-8892-638c4ff40052" />
whisker plot indication a good sequencing .
<img width="1008" height="743" alt="Screenshot 2025-10-31 161224" src="https://github.com/user-attachments/assets/2f7c141f-26cb-488a-aa2a-eb5193034912" />
<img width="984" height="735" alt="Screenshot 2025-10-31 161236" src="https://github.com/user-attachments/assets/e03c7068-26be-487f-8556-760f8fd0018a" />
usually the high peak of specific base pairs seen in the starting is a result of primer sequences attached, adapters and amplification for library preparation, a common glitch which is acceptable in sequencing analyses and can be trimmed using trimmonatic.
<img width="979" height="744" alt="Screenshot 2025-10-31 161246" src="https://github.com/user-attachments/assets/2eaa5c58-aba1-4dd6-995f-dc336625f26e" />
<img width="979" height="747" alt="Screenshot 2025-10-31 161257" src="https://github.com/user-attachments/assets/3e911f46-d4b4-4c8f-9f63-1a9d0a98aa68" />
here we can see that all the sequences are indentified as one of the four base pairs, and hence no N content.
<img width="962" height="209" alt="Screenshot 2025-10-31 161354" src="https://github.com/user-attachments/assets/d9996659-e75a-448d-b033-ee3c61fdf1e8" />
the sequence duplicates are usually a result of amplification process and priming for library preparation.
```bash
sudo apt install multiqc
```
```bash
multiqc fastqc_results/ -o multiqc_report/
```
the multiqc reports are as follows:
<img width="1568" height="636" alt="image" src="https://github.com/user-attachments/assets/c9f71e84-e32c-4145-9664-7b0a2c9a55a3" />
image showing us a representation of multiqc reports and read counts per sample, a great tool to comparatively analyse the samples with each other simultaniously.
# Trimming using Trimmonatic and Quality check
usually trimmimg is required to removed any adapter sequences and repeated reads . but as we saw in our case, the data was already clean and ready to be used ahead. This can be further expressed by comparing the fastqc results before and after trimming.
```bash
java -jar trimmomatic-0.40.jar SE -threads 4 fastq/LNCAP_Hypoxia_S1.fastq.gz fastq/LNCAP_Hypoxia_S1_trimmed.fastq TRAILING:10 -phred33
```
```bash
fastqc fastq/*.fastq.gz -o fastqc_results/ --threads
```
# Downloading human reference genome
```bash
wget ftp://ftp.ensembl.org/pub/release-114/gtf/homo_sapiens/Homo_sapiens.GRCh38.114.gtf.gz
```
unzipping the file
```bash
gunzip Homo_sapiens.GRCh38.114.gtf.gz
```
# Alignment using Hisat2 and Samtools
Usually, Hisat2 and STAR are the softwares used for alignment of reads to the reference genome . In our case we have used Hisat2 as STAR is a larger application to download and might even take days to process.
```bash
sudo apt install hisat2
```
```bash
sudo apt install samtools
```
the step coming forth is an automated pipeline used to call upon each sample file and align them over human genome and producing the .bam files. 
```bash
./hisat2alignment.sh
```
<img width="1088" height="303" alt="Screenshot 2025-10-03 235036" src="https://github.com/user-attachments/assets/89dacd69-cc39-4a0b-8d79-99b80c2a0f78" />
# Installing Qualimap
```bash
conda install -c bioconda qualimap
```
```bash
./qualimap_v2.3qualimap rnaseq -bam yourpath/sample.bam -gtf Homo_sapiens.GRCh38.114.gtf –outdir yourpath/outputfolder –java-mem-size=8G
```
# Generating freature counts using subread
```bash
sudo apt install subread
```
```bash
mkdir -p quants
```
```bash
./featurecounts.sh
```
The above code runs a script to generate feature counts for each file in a loop fashion.
<img width="1156" height="920" alt="Screenshot 2025-10-04 233323" src="https://github.com/user-attachments/assets/04afd812-7e18-4dcd-a5f1-303f2d4bebb0" />

# Generating counts matrix from feature counts 
This entire section was done on google colab. this file will be links in the codes section.

# Gene annotation of data
assigning a gene prototype is an important step before R analysis.

## R analysis for differential gene expression. 
The entire R analysis is pinned under the R section. Results and command line are available there.






