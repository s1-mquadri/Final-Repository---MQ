# Final-Repository---MQ
This is a repository that follows the steps taken to create the results discussed in my research paper about the NOC4L gene from the CBF/Mak21 family. 

To begin, you need to make sure your instance is running on EC2 and to log in via ssh. 


Then, create a new working directory. 
```
mkdir ~/globins2
```

Now go to this folder: 
```
cd ~/globins2
```

Next, download your query protein. For me, it was “NP_076983.1”
```
ncbi-acc-download -F fasta -m protein "NP_076983.1"
```
Then, perform the BLAST search. 
```
blastp -db ../allprotein.fas -query NP_076983.1.fa -outfmt 0 -max_hsps 1 -out globins2.blastp.typical.out
seqkit grep --pattern-file ~/globins2/globins.blastp.detail.filtered.out ~/globins2/NP_076983.1.fa | seqkit grep -v -p "carpio" > ~/globins2/globins2.homologs.fas
```

To look through the results, we can use “less.”
```
less globins2.blastp.typical.out
```
Now, we will create an easier to read output for the same BLAST output. We are adding the “outfmt” flag that specifies a specific output formate 
```
blastp -db ../allprotein.fas -query NP_076983.1.fa  -outfmt "6 sseqid pident length mismatch gapopen evalue bitscore pident stitle"  -max_hsps 1 -out globins2.blastp.detail.out
```
Now, let’s filter the BLAST output for high-scoring putative homologs. We’ll require the e-value to be less than 1e-30. To do this, we will use the awk program, which will select particular records in a file and perform operations upon them.

```
awk '{if ($6< 1e-30)print $1 }' globins2.blastp.detail.out > globins2.blastp.detail.filtered.out
```
To count the total number of hits in the BLAST results after the filter, we will use the “wc command.” 
```
wc -l globins2.blastp.detail.filtered.out
```

To find out how many paralogs are in each species, we will use this command. 
```
grep -o -E "^[A-Z]\.[a-z]+" globins2.blastp.detail.filtered.out  | sort | uniq -c
```
To perform the global multiple sequence alignment, using MUSCLE, we will use this code. MUSCLE is a common alignment program. Aligning all the sequences will show which positions are homologous to each other, the insertions/deletions, etc. 
```
muscle -align ~/globins2/globins2.homologs.fas -output ~/globins2/globins2.homologs.al.fas
```
Now, use the R-package msa and a script to create a pdf of your alignment. 
```
 Rscript --vanilla ~/plotMSA.R  ~/globins2/globins2.homologs.al.fas
```
To calculate the length of the alignment (with column gaps), use this code. 
```
alignbuddy  -al  ~/globins2/globins2.homologs.al.fas
```
To calculate the length (without column gaps), use this code. 
```
alignbuddy -trm all  ~/globins2/globins2.homologs.al.fas | alignbuddy  -al
```
To calculate the average percent identity, we will use two methods, t_coffee, and align buddy. First, t_coffee. 
```
t_coffee -other_pg seq_reformat -in ~/globins2/globins2.homologs.al.fas -output sim
```
Now, align buddy. 
```
 alignbuddy -pi ~/globins2/globins2.homologs.al.fas | awk ' (NR>2)  { for (i=2;i<=NF  ;i++){ sum+=$i;num++} }
     END{ print(100*sum/num) } '
```
We will use IQ-TREE, a software to predict optimal phyloggenetic tree based on our sequence alignment. This software will first calculate the optimal amino acid substitution model and amino acid frequencies, then it will perform a tree search and will estimate branch lengths. The first line will remove any duplicate label tags, and the second line will run the IQ-TREE estimation. 
```
sed 's/ /_/g'  ~/globins2/globins2.homologs.al.fas | seqkit grep -v -r -p "dupelabel" >  ~/globins2/globins2.homologsf.al.fas
```
```
iqtree -s ~/globins2/globins2.homologsf.al.fas -bb 1000 -nt 2 
```
To look at the unrooted tree we will use the following command. 
```
nw_display ~/globins2/globins2.homologsf.al.fas.treefile
```
To look at this same output in a graphic display, use this command. 
```
Rscript --vanilla ~/plotUnrooted.R  ~/globins2/globins2.homologsf.al.fas.treefile ~/globins2/globins2.homologsf.al.fas.treefile.pdf 0.4 15
```
Now, we will use midpoint rooting to root the tree so that we can understand the oldest divergence event. 
```
gotree reroot midpoint -i ~/globins2/globins2.homologsf.al.fas.treefile -o ~/globins2/globins2.homologsf.al.mid.treefile
```
To view this tree as a pdf, use these commands. 
```
nw_order -c n ~/globins2/globins2.homologsf.al.mid.treefile  | nw_display -
```
```
nw_order -c n ~/globins2/globins2.homologsf.al.mid.treefile | nw_display -w 1000 -b 'opacity:0' -s  >  ~/globins2/globins2.homologsf.al.mid.treefile.svg -
```
```
convert  ~/globins2/globins2.homologsf.al.mid.treefile.svg  ~/globins2/globins2.homologsf.al.mid.treefile.pdf
```
Now we will reconcile the gene tree and species tree. 
```
java -jar ~/tools/Notung-3.0_24-beta/Notung-3.0_24-beta.jar -s ~/species.tre -g ~/globins2/globins2.homologsf.al.mid.treefile --reconcile --speciestag prefix --savepng --events --outputdir ~/globins2/
```
To look at this output, use the following command. 
```
less globins2.homologsf.al.mid.treefile.rec.events.txt

```
To generate a RecPhyloXML object and to view the gene-tree species via thirdkind, use this command. 
```
python2.7 ~/tools/recPhyloXML/python/NOTUNGtoRecPhyloXML.py -g ~/globins/2globins2.homologsf.al. mid.treefile.rec.ntg --include.species
```
To create the reconciliation graphic, use this command. 
```
thirdkind -Iie -D 40 -f ~/globins2/globins2.homologsf.al.mid.treefile.rec.ntg.xml -o  ~/globins2/globins2.homologsf.al.mid.treefile.rec.svg
```
Convert this to a pdf. 
```
convert  -density 150 ~/globins2/globins2.homologsf.al.mid.treefile.rec.svg ~/globins2/globins2.homologsf.pruned.treefile.rec.pdf
```
Now, let’s run the RPS-BLAST to identify Pfam domains in our protein sequences. 
```
rpsblast -query ~/globins2/globins2.homologs.fas -db ~/data/Pfam/Pfam -out ~/globins2/globins2.rps-blast.out  -outfmt "6 qseqid qlen qstart qend evalue stitle" -evalue 0.5
```
Now, we will run the Rscript to plot the predicted Pfam domains on the phylogeny. 
```
Rscript  --vanilla ~/lab08-$MYGIT/plotTreeAndDomains.r ~/globins2/globins2.homologsf.al.mid.treefile ~/globins2/globins2.rps-blast.out ~/globins2/globins2.tree.rps.pdf
```



