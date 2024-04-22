# Get TajimaD with 10K window
vcftools --vcf POP1.vcf --TajimaD 10000 --out POP1.10K.TajimaD
# Get per window > 10 SNPs 
awk '$3 >= 10' POP1.10k.Tajima.D > filtered_POP1.10k.Tajima.D
# Annotated with BEDTools to get bottom 5%
sed -e '/nan/d' filtered_POP1.10k.Tajima.D > TajimaD.1.txt
sort -k4gr  TajimaD.1.txt | tail -n $(($(wc -l <  TajimaD.1.txt) *5 / 100)) > TajimaD.2.txt
awk '{print $1, $2, $2 + 10000}' TajimaD.2.txt > TajimaD.3.txt 
sed -i '$ d' TajimaD.3.txt
sed 's/ /\t/g'  TajimaD.3.txt >  TajimaD.4.txt
sort -k1,1n -k2,2n   TajimaD.4.txt >   TajimaD.5.txt
bedtools intersect -a Reference.Genome.sort.gff3 -b   TajimaD.5.txt  -wa >   TajimaD.6.txt 
awk -F' ' '{print $9}'  TajimaD.6.txt | grep -oP '(?<=Parent=)[^,]*' | sed  's/;//g' >  TajimaD.7.txt
sort   TajimaD.7.txt | uniq >   TajimaD.8.txt
grep -v 'Name'  TajimaD.8.txt | grep '\.t1$' > TajimaD.bottom5%.final.GeneID.uniq.txt
