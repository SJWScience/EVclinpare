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
awk 'FNR==1{print ""}1' *.fasta > catfile_allgenomes.fasta

makeblastdb -in catfile_allgenomes.fasta -dbtype nucl -title Specific_Genomes -out Genome_database

```

