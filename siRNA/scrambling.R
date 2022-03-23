
library(seqRFLP)
sequ <- "ACAAGGTGAAGACGCTCAA"
separate.seq.letters <- unlist(strsplit(sequ, ""))
l=replicate(100000, sample(1:nchar(sequ)), simplify=F)
n=lapply(l,function(x){paste(separate.seq.letters[x], collapse = "")})
m=as.data.frame(unique(unlist(n)))
write.table(dataframe2fas(m, file = './scr_multi.fa'))

system('makeblastdb -in ./Homo_sapiens.GRCh38.cdna.all.fa -parse_seqids -out ./Homo_sapiens_GRCh38_cdna -dbtype nucl')
system('blastn \
-num_threads 8 \
-task blastn-short \
-max_target_seqs 1 \
-evalue 10000 \
-query ./scr_multi.fa \
-db ./Homo_sapiens_GRCh38_cdna \
-out ./scr_blsast.tab \
-outfmt 6 ')

jd=read.table('./scr_blsast.tab')
max(jd$V11)
length(unique(jd$V1))
jdt=tapply(jd$V11,jd$V1,min)
jdt=as.data.frame(jdt)
jdt=cbind(rownames(jdt),jdt)
jdm=jdt[jdt$jdt==283,]
colnames(jdm)=c('V1','V11')
m=merge(jdm,jd,by=c('V1','V11'))
write.table(m,'./non_homolog_scr_blsast.tab',quote=F,sep='\t')



