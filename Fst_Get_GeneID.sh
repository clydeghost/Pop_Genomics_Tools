# Get Fst with 10k window or  Get Fst with POPgenWindows,different seps with two methods
vcftools --vcf All.vcf --weir-fst-pop POP1.txt --weir-fst-pop POP2.txt --fst-window-size 10000 --out POP1_POP2_Fst_10k
# Get per windows > 10 SNPs
awk -F ',' '$5 >= 10' POP1_POP2.10k.5k.Fst.csv  > filtered_POP1_POP2.10k.5k.Fst.csv 
# Annotated with BEDTools
sort -t ',' -k9gr filtered_POP1_POP2.10k.5k.Fst.csv  | head -n $(($(wc -l < POP1_POP2.10k.5k.Fst.csv) * 5 / 100)) > Fst.1.txt
awk -F ',' '{print$1,$2,$3}' Fst.1.txt  > Fst.2.txt 
sed 's/ /\t/g' Fst.2.txt > Fst.3.txt
sort -k1,1n -k2,2n  Fst.3.txt >  Fst.4.txt
bedtools intersect -a Reference.Genome.sort.gff3 -b  Fst.4.txt  -wa >  Fst.5.txt 
awk -F' ' '{print $9}' Fst.5.txt | grep -oP '(?<=Parent=)[^,]*' | sed  's/;//g' > Fst.6.txt
sort  Fst.6.txt | uniq >  Fst.7.txt
grep -v 'Name' Fst.7.txt | grep '\.t1$' > Fst.top5%.final.GeneID.uniq.txt
