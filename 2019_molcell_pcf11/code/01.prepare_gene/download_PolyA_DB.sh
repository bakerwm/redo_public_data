
## download PolyA_DB human
wget --no-check-certificate https://exon.apps.wistar.org/polya_db/v3/download/3.2/human_pas.zip
unzip human_pas.zip # output: human.PAS.txt
awk -F"\t" 'OFS="\t"{print $2,$3-1,$3,$1,"254",$4}' human.PAS.txt | sed '1d' > human_pas.hg19.bed
awk -F"\t" 'OFS="\t"{print $2,$3-1,$3,$1,"254",$4,$14}' human.PAS.txt | grep -P "3'UTR\(L\)|3'UTR\(S\)" | cut -f1-6 > human_pas.hg19.3UTR.bed
wc -l human*bed
# expect
#    16781 human_pas.hg19.3UTR.bed
#   311594 human_pas.hg19.bed

