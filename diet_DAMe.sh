#Running a modified DAMe pipeline from https://github.com/shyamsg/DAMe 
#We run each marker in their seperate directory

#First create a sub-directory for each sample/library and keep the raw read in these. 
mkdir pool1 
mkdir pool2 

#Now enter the the sub-directory and run the following commands for each of the two pools.
cd pool1
#We use AdapterRemoval to trim, remove adaptors and merge paired reads
module load AdapterRemoval

xsbatch -c 1 -- AdapterRemoval --file1 ${sample}_1_R1.fastq.gz --file2 ${sample}_1_R2.fastq.gz --mm 0.05 --minlength 100 --shift 5 --basename pool1/pool1_merged --trimns --trimqualities --qualitybase 33 --minquality 28 --minalignmentlength 50 --collapse

#Combine the merged collapsed files
cat pool1_merged.collapsed pool1_merged.collapsed.truncated > pool1_merged.fastq

#Move to project folder and the we'll sort the reads by primer and tag.
module load python/v2.7.12
module load DAMe

xsbatch -c 1 -- python DAMe/v0.9/bin/sort.py -fq /pool1/pool1_merged.fastq -p ${sample}/${sample}_primers.txt -t ${sample}/${sample}_tags.txt

 #We can get a summary of the tag combinations within each pool with the following command:
xsbatch -c 1 -- python /DAMe/v0.9/bin/splitSummaryByPSInfo.py -p ${sample}/${sample}_PSinfo.txt -l 1 -s /pool1/SummaryCounts.txt -o /pool1/SummaryCounts_split_pool1.txt

#Now repeat the above commands for pool2 
#Now we want to go back to the project folder and sort the reads into their respective biological samples

xsbatch -c 1 -- python /DAMe/v0.9/bin/filter.py -psInfo /${sample}_PSinfo.txt -x 3 -y 2 -p 2 -t 2 -l 90
#The -l parameter needs to be set according to marker length, we used the following Coleop -l 90, Mam16S -l 80, Trac01 -l 250, Vert01 -l 80 and Zeale -l 140

#Now we need to convert the filtered reads into a format we can use in the next step
xsbatch -c 1 -- python /DAMe/v0.9/bin/convertToUSearch.py -i /FilteredReads.fna

#Then we will create OTU clusters using sumaclust
module load sumaclust/v1.0.20

xsbatch -c 1 -- sumaclust -e FilteredReads.forsumaclust.fna -F OTUs_97_sumaclust.fna    #We use default settings and get a cut off at 97%

#Normalise sequence numbers and create an OTU table
xsbatch -c 1 -- python /DAMe/v0.9/bin/tabulateSumaclust.py -i OTUs_97_sumaclust.fna -s 50000 -o table_${sample}_97.txt -blast

#Assign Taxonomy with blast
module load blast+/v2.6.0

xsbatch -c 1 --mem-per-cpu 35000 -- blastn -query table_${sample}_97.txt.blast.txt -out blast_${sample}_97.output.txt -db nt

#Export the table_${sample}_97.txt and blast_${sample}_97.output.txt for data analysis. 