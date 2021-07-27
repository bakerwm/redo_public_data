# Download PolyASite Human PAS annotation to PWD

wget https://polyasite.unibas.ch/download/atlas/2.0/GRCh38.96/atlas.clusters.2.0.GRCh38.96.bed.gz
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz
# hg38
zcat atlas.clusters.2.0.GRCh38.96.bed | cut -f1-6 | sed 's/^/chr/' > pas.2.0.hg38.bed
zcat atlas.clusters.2.0.GRCh38.96.bed | awk '$10=="TE"{print $0}' | cut -f1-6 | sed 's/^/chr/' > pas.2.0.hg38.TE.bed
# hg19
CrossMap.py bed hg38ToHg19.over.chain.gz pas.2.0.hg38.bed pas.2.0.hg19.bed
CrossMap.py bed hg38ToHg19.over.chain.gz pas.2.0.hg38.TE.bed pas.2.0.hg19.TE.bed
# output
wc -l pas.2.0*bed 
# expect
#   143588 pas.2.0.hg19.TE.bed
#   567729 pas.2.0.hg19.bed
#   569005 pas.2.0.hg38.bed

