library(gplots)

stem_score=read.table('stem_score.txt',row.names=1,header=F,sep='\t')
ac_score=read.table('AC_score.txt',row.names=1,header=F,sep='\t')
oc_score=read.table('OC_score.txt',row.names=1,header=F,sep='\t')
IDENT=read.table('IDENT.txt',row.names=1,header=F,sep='\t')
B= which(rownames(stem_score) %in% rownames(IDENT))
stem_score_b=stem_score[B,]
ac_score_b=ac_score[B,]
oc_score_b=oc_score[B,]


pdf('STEM_AC_OC.pdf',width=10,height=10)
COMBINE=c()
##########STEM###########
stem_score_b=stem_score[B]
med_stem_score=c()
wilcox_test_p=c()
i=0
while(i<length(table(IDENT))){
med_stem_score=c(med_stem_score, median(stem_score_b[which(IDENT==i)]))
wilcox_test_p=c(wilcox_test_p,wilcox.test(stem_score_b[which(IDENT==i)],stem_score_b)$p.value ) 
i=i+1
}
neg_log_2_p=-log(wilcox_test_p,2)
neg_log_2_adjp= -log(p.adjust(wilcox_test_p,method='fdr'),2)
direction= (med_stem_score - median(stem_score_b))/abs((med_stem_score - median(stem_score_b)))
#plot(neg_log_2_p*direction,pch=16)
#abline(h=0)
#abline(h=-log(0.05,2),lty=3)
#abline(h=log(0.05,2),lty=3)
COL=rep('black',length(neg_log_2_adjp))
COL[which(neg_log_2_adjp*direction > -log(0.05,2))]='red'
plot(x=c(0: (length(table(IDENT))-1) ),y=neg_log_2_adjp*direction,pch=16,ylim=c(-150,50),col=COL,cex=3,main='STEM',xlab='cluster number',ylab='fdr')
abline(h=0)
abline(h=-log(0.05,2),lty=3)
abline(h=log(0.05,2),lty=3)

COMBINE=cbind(COMBINE,neg_log_2_adjp*direction)
#################################

##########AC###########
stem_score_b=ac_score_b
med_stem_score=c()
wilcox_test_p=c()
i=0
while(i<length(table(IDENT))){
med_stem_score=c(med_stem_score, median(stem_score_b[which(IDENT==i)]))
wilcox_test_p=c(wilcox_test_p,wilcox.test(stem_score_b[which(IDENT==i)],stem_score_b)$p.value ) 
i=i+1
}
neg_log_2_p=-log(wilcox_test_p,2)
neg_log_2_adjp= -log(p.adjust(wilcox_test_p,method='fdr'),2)
direction= (med_stem_score - median(stem_score_b))/abs((med_stem_score - median(stem_score_b)))
#plot(neg_log_2_p*direction,pch=16)
#abline(h=0)
#abline(h=-log(0.05,2),lty=3)
#abline(h=log(0.05,2),lty=3)
COL=rep('black',length(neg_log_2_adjp))
COL[which(neg_log_2_adjp*direction > -log(0.05,2))]='red'
plot(x=c(0: (length(table(IDENT))-1) ),y=neg_log_2_adjp*direction,pch=16,col=COL,cex=3,main='AC',xlab='cluster number',ylab='fdr')
abline(h=0)
abline(h=-log(0.05,2),lty=3)
abline(h=log(0.05,2),lty=3)
COMBINE=cbind(COMBINE,neg_log_2_adjp*direction)
#################################

##########OC###########
stem_score_b=oc_score_b
med_stem_score=c()
wilcox_test_p=c()
i=0
while(i<length(table(IDENT))){
med_stem_score=c(med_stem_score, median(stem_score_b[which(IDENT==i)]))
wilcox_test_p=c(wilcox_test_p,wilcox.test(stem_score_b[which(IDENT==i)],stem_score_b)$p.value ) 
i=i+1
}
neg_log_2_p=-log(wilcox_test_p,2)
neg_log_2_adjp= -log(p.adjust(wilcox_test_p,method='fdr'),2)
direction= (med_stem_score - median(stem_score_b))/abs((med_stem_score - median(stem_score_b)))
#plot(neg_log_2_p*direction,pch=16)
#abline(h=0)
#abline(h=-log(0.05,2),lty=3)
#abline(h=log(0.05,2),lty=3)
COL=rep('black',length(neg_log_2_adjp))
COL[which(neg_log_2_adjp*direction > -log(0.05,2))]='red'
plot(x=c(0: (length(table(IDENT))-1) ),y=neg_log_2_adjp*direction,pch=16,col=COL,cex=3,main='OC',xlab='cluster number',ylab='fdr')
abline(h=0)
abline(h=-log(0.05,2),lty=3)
abline(h=log(0.05,2),lty=3)
COMBINE=cbind(COMBINE,neg_log_2_adjp*direction)
#################################
dev.off()



#################################



