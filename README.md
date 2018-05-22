# EVclinpare
Basic outline for comparison of discovered mutations in lab evolution experiments to environmental and clinical isolates of Pseudomonas aeruginosa


## Creating databases for searching

> First step for these procedure is to create a database(s) of clinical or environmental isolates. Typically these are not well assembled or are in a state of raw contigs output from an assembly program. So usually some formatting of these files is recommended for the cleanest results. This can be achieved by using two small scripts that will rename the file headers from generic contig headers and will remove any internal fasta headers as they can interfere with downstream analyses.

```bash
for i
in $(ls *.fa*); do grep -i '^[ACGTURYKMSWBDHVNX\-]*$' "$i" | tr -d '\n' | fold -w 60 > /"$i".fmthead
done

```

> and

```bash
for i 
in $(ls *.fmthead); 
do
sed '1 i\>'"$i"'' "$i" > "$i".fasta
done
```

> This is pretty dirty code but it does the job of formatting these files. There are probably much simplier and quicker ways of doing this but currently this is how i am doing it.

> Next step is to create the blast databases for each file type (if using env and clinical genomes) or just to create one database if only using one type.

```bash
awk 'FNR==1{print ""}1' *.fasta > catfile_allgenomes.forDB

makeblastdb -in catfile_allgenomes.fasta -dbtype nucl -title Specific_Genomes -out Genome_database

```
## Searching the database

> For a more hands on manual analysis (useful for a handful of genes) then you can do this analysis gene by gene - shown below. NOTE THIS IS SPECIFIC TO BIOC DEPT at UoO - other programs etc will not be around elsewhere

```bash
export PATH=/Volumes/userdata/staff_groups/lamontlab/Documents/BLAST/output_files/mview-1.60.1/bin/:$PATH

tblastn -query input_genes/geneA.in -db Genome_database -out output_results/geneA.blast -num_threads 100 -max_hsps 1 -num_alignments 400

mview -top 500 -in blast output_results/geneA.blast -out fasta > output_results/geneA.blast.fasta

snp-sites -v output_results/geneA.blast.fasta -o output_results/geneA_output_variants.vcf

```

> If you have lots of genes of interest then you can use a script to parse them all from one file to another this is EVvlinpare

```bash
for i
in $(ls *.in | rev | cut -c 4- | rev)
do
tblastn -query "$i".in -db Genome_database -out output_results/"$i".blast -num_threads 100 -max_hsps 1 -num_alignments 400
echo "$i" MVIEWS
mview -top 500 -in blast output_results/"$i".blast -out fasta > output_results/"$i".blast.fasta
#echo "$i" SNP-SITES
snp-sites -v output_results/"$i".blast.fasta -o output_results/"$i".output_variants.vcf
done

```

>Some genes may not work with this script - this in my experience is due to there being very large amounts of difference between the reference sequence and the database genomes being used. However this is very rare in my experience. To get around this if it does happen there is some slight alterations you can make to the code.

```bash
tblastn -query geneA.in -db Genome_database -out geneA.blastout -num_threads 100 -max_hsps 1 -num_alignments 400 -outfmt '6 qseqid sseqid sseq'

awk 'BEGIN { OFS = "\n" } { print ">"$2, $3 }' geneA.blastout | cat > geneA.fasta_for_aln

clustalo -i geneA.fasta_for_aln -o geneA.omegaALN -v --threads=100

snp-sites -v geneA.omegaALN -o geneA.output_variants.vcf

```

## Analysing results

> The output from the above pipeline will give you a .vcf (variant call format) file. These can be opened by any text editor or can be shown in excel (recommended for inexperienced users). Where it will show the postion within the gene that has been mutated along with the reference residue and variant residue. There will also be a column for each genome you gave it with an integer corresponding to presence/absence of mutation - 1 being mutation corresponding to the first listed variant of that row, 2 is 2nd mutation in the same row . . . and so on. Most of the time there will only be one variant per postion, however, there are genes that are quite variable or have specific residues that are commonly mutated. In cases like this you will see multiple variants for the same potion and a genome may have a number above 1. If there is a * or an X it means there was either no sequence called at this postion (possible deletion) or unable to confidently call a residue (usually from the assembly being unsure about a specific base and assigning an ambiguous base to it. If these are issues for your analysis go back to the raw blast output and open it in a text editior and examine the genome(s) of interest to see if there is a deletion or if something else is happening. In my experience its sorta a coinflip as to if its a deletion or an assembly thing. Quite often it was the beginning or the end of a contig and thus the blast missed this region because the genomes contigs are not necessarily in order.

> Input variants can be parsed into a program to call putative function of mutations - EG: provean. To get this working on command line is a bit of a nightmare in my experience, but once its done it is actually quite good.

**Playing around with automated parsing of mutations in command line to streamline the entire process

```bash
 awk '/#CHROM/ ,EOF { print $4,$2,$5 }' fusA1.CFALTERNATE_ONLYCF.vcf > SAMTRIAL.txt

```
> This is pretty cool it prints just the variants in the right order, now i just need a way to handle muliple entries that have the commas etc.

```bash
awk '{ n=split($3,arr,","); for(i=1;i<=n;i++) print $1,$2,arr[i] }' SAMTRIAL.txt | tr -d "[:blank:]" 

awk '{ n=split($3,arr,","); for(i=1;i<=n;i++) print $1,$2,arr[i] }' SAMTRIAL.txt | tr -d "[:blank:]" | sed '1d' > OUTPUT.var

```

**How to automatically extract variants from the VCF files

> After some work i finally got it going

```bash
awk '/#CHROM/ ,EOF { print $4,$2,$5 }' geneA.vcf | awk '{ n=split($3,arr,","); for(i=1;i<=n;i++) print $1,$2,arr[i] }' | tr -d "[:blank:]" | sed '1d' > OUTPUT_geneA.var

```
