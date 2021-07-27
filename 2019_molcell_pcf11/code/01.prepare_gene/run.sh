

# download data
bash download_PolyASite.sh 
bash download_PolyA_DB.sh 
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz

# extract info
python filt_gtf.py -i gencode.v19.annotation.gtf.gz -o to_PolyASite_TE -b pas.2.0.hg19.TE.bed -f gene -x --size 5000 --flank 6000
python filt_gtf.py -i gencode.v19.annotation.gtf.gz -o to_PolyASite_full -b pas.2.0.hg19.bed -f gene -x --size 5000 --flank 6000
python filt_gtf.py -i gencode.v19.annotation.gtf.gz -o to_PolyA_DB_3UTR -b human_pas.hg19.3UTR.bed -f gene -x --size 5000 --flank 6000 
python filt_gtf.py -i gencode.v19.annotation.gtf.gz -o to_PolyA_DB_full -b human_pas.hg19.bed -f gene -x --size 5000 --flank 6000 
