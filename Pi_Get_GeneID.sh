# Get Pi with Vcftools
vcftools --vcf ALL.vcf --keep POP1.txt --window-pi 10000 --out POP1.10k.windowed.pi
# Get per windows > 10 SNPs
awk '$4 >= 10' POP1.10k.windowed.pi > filtered_POP1.10k.windowed.pi
awk '$4 >= 10' POP2.10k.windowed.pi > filtered_POP2.10k.windowed.pi
# Get same windows
awk 'FNR==NR{a[$1,$2];next}($1,$2) in a' filtered_POP1.10k.windowed.pi   filtered_POP1.10k.windowed.pi > POP1.txt
awk 'FNR==NR{a[$1,$2];next}($1,$2) in a' filtered_POP2.10k.windowed.pi   filtered_POP2.10k.windowed.pi > POP2.txt
# Fet Pi Radio
awk 'FNR==NR{a[FNR]=$5;next}(FNR>1){$5=a[FNR]/$5}1' POP1.txt  POP2.txt > POP1_POP2.10k.PiRadio.txt
# Annotated with BEDTools
sort -k5gr  POP1_POP2.10k.PiRadio.txt | head -n $(($(wc -l <  POP1_POP2.10k.PiRadio.txt) *5 / 100)) > Pi.1.txt
awk '{print $1, $2, $3}' Pi.1.txt > Pi.2.txt
sed 's/ /\t/g'  Pi.2.txt > Pi.3.txt
sort -k1,1n -k2,2n   Pi.3.txt > Pi.4.txt
bedtools intersect -a Reference.Genome.sort.gff3  -b   Pi.4.txt  -wa >  Pi.5.txt
awk -F' ' '{print $9}'  Pi.5.txt | grep -oP '(?<=Parent=)[^,]*' | sed  's/;//g' >  Pi.6.txt
sort   Pi.6.txt | uniq >  Pi.7.txt
grep -v 'Name'  Pi.7.txt | grep '\.t1$' > Pi.top5%.final.GeneID.uniq.txt
