# Get biaallelic vcf
bcftools view --types snps -m 2 -M 2 input.vcf -o  ALL.biallelic.vcf
# Run XPEHH 
for i in $(cat ../Chr.txt);
do
vcftools --gzvcf ALL.biallelic.vcf   --max-missing 1  --recode --recode-INFO-all --chr chr${i} --out ALL.chr${i}
vcftools --vcf ALL.chr${i}.recode.vcf  --plink --out ALL.chr${i};

awk 'BEGIN{OFS=" "} {print 1,".",$4/1000000,$4}' ALL.chr${i}.map > ALL.chr${i}.map.distance;

vcftools --vcf ALL.chr${i}.recode.vcf  --recode --recode-INFO-all --chr chr${i} --keep POP1id.txt  --out ALL.chr${i}.SJLS;
vcftools --vcf ALL.chr${i}.recode.vcf  --recode --recode-INFO-all --chr chr${i} --keep POP2id.txt --out ALL.chr${i}.KX;

selscan-2.0.0/src/selscan --xpehh --vcf ALL.chr${i}.POP1.recode.vcf --vcf-ref ALL.chr${i}.POP2.recode.vcf --map ALL.chr${i}.map.distance --out POP1_POP2.chr${i}  --threads 24 --ehh-win 10000;
wait
awk  '{print '${i}',$2,$3,$4,$5,$6,$7,$8}' POP1_POP2.chr${i}.xpehh.out > POP1_POP2.chr${i}.xpehh.csv;
sed -i 's/ /\t/g' POP1_POP2.chr${i}.xpehh.csv;
selscan-2.0.0/src/norm --xpehh --files POP1_POP2.chr${i}.xpehh.csv --bp-win --winsize 10000
done
# Get Reference.Genome.sort.gff3
bedtools sort Reference.Genome.gff3 > Reference.Genome.sort.gff3
# Annotated through BEDTools to get top 5%
sort -k9gr  POP1_POP2.xpehh.csv.norm.10kb.windows | head -n $(($(wc -l < POP1_POP2.xpehh.csv.norm.10kb.windows) * 5 / 100)) > XPEHH.top5%.txt
awk '{print$1,$2,$3}'  XPEHH.top5%.txt >  XPEHH.top5%.1.txt
awk '{$1 = "Fp"$1}1' XPEHH.top5%.1.txt > XPEHH.top5%.2.txt
sed 's/ /\t/g' XPEHH.top5%.2.txt > XPEHH.top5%.3.txt
sort -k1,1n -k2,2n  XPEHH.top5%.3.txt >  XPEHH.top5%.3.sort.txt
bedtools intersect -a Reference.Genome.sort.gff3 -b  XPEHH.top5%.3.sort.txt  -wa >  XPEHH.top5%.annotated.txt
awk -F ' ' '{print $9}' XPEHH.top5%.annotated.txt | grep -oP '(?<=Parent=)[^,]*' | sed  's/;//g' > XPEHH.top5%.anno.GeneID.txt
sort  XPEHH.top5%.anno.GeneID.txt | uniq >  XPEHH.top5%.anno.GeneID.uniq.txt
grep -v 'Name' XPEHH.top5%.anno.GeneID.uniq.txt | grep '\.t1$' > XPEHH.top5%.final.GeneID.uniq.txt