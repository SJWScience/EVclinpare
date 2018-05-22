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










