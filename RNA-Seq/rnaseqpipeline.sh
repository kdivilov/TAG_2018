#RNA-Seq fastq files: https://www.ncbi.nlm.nih.gov/bioproject/281110

#HISAT2: https://ccb.jhu.edu/software/hisat2/index.shtml (Note that a newer version has been released and the code below should be changed to reflect this if using the newer version)
#StringTie: https://ccb.jhu.edu/software/stringtie/
#Ballgown: https://github.com/alyssafrazee/ballgown

#Genome assembly file: https://urgi.versailles.inra.fr/Species/Vitis/Data-Sequences/Genome-sequences
#Annotation file: https://urgi.versailles.inra.fr/Species/Vitis/Annotations

#Create index
hisat2-2.0.5/hisat2-build -p 8 ref/12Xv2_grapevine_genome_assembly.fa ref/Vitisvinifera12X2

#Align reads
for filename in *.fastq.gz; do
    echo "$filename"
    hisat2-2.0.5/hisat2 -p 8 --dta -x ref/Vitisvinifera12X2 -U "$filename" -S "hisatoutput/$(basename #"$filename" .fastq.gz)_hisat.sam" 
done

#SAM->BAMsorted
for filename in hisatoutput/*.sam; do
    echo "$filename"
    samtools sort -@ 8 -o "hisatoutput/$(basename "$filename" .sam).bam" "$filename"
done

#StringTie

for filename in hisatoutput/*.bam; do
    echo "$filename"
    stringtie-1.3.3b.Linux_x86_64/stringtie "$filename" -o "stringtieoutput/$(basename "$filename" .bam).gtf" -p 8 -G ref/all_V1_NR_mrna_new_posit_on_assembly_v2_without_chrUkn.gff3
done

stringtie-1.3.3b.Linux_x86_64/stringtie --merge -p 8 -G ref/all_V1_NR_mrna_new_posit_on_assembly_v2_without_chrUkn.gff3 -o  stringtieoutput/allmerged.gtf mergelist.txt

#Setup for Ballgown analysis

for filename in stringtieoutput/*.gtf; do
    echo "$filename"
    stringtie-1.3.3b.Linux_x86_64/stringtie "hisatoutput/$(basename "$filename" .gtf).bam" -e -B -p 8 -G stringtieoutput/allmerged.gtf -o "ballgown/BG_$(basename "$filename" .gtf)/BG_$(basename "$filename" .gtf).gtf" 
done



