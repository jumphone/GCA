library(gplots)

exp_data=as.matrix(read.table('GSE70630_OG_processed_data_v2.txt.cleaned.txt',header=T,row.names=1))
tf_ident=read.table('IDENT.txt')


EXP = CreateSeuratObject(raw.data = b, min.cells = 0, min.genes=0)

EXP <- AddMetaData(object = EXP, metadata = stem_score_b, col.name = "stem.score")
EXP@meta.data$stem.score=stem_score_b

EXP=NormalizeData(object = EXP, normalization.method = "LogNormalize", scale.factor = 10000)

EXP = ScaleData(object = EXP,, genes.use = all_gene)

PCNUM=40
EXP <- RunPCA(object = EXP, pc.genes = all_gene, do.print = TRUE, pcs.print = 1:5,    genes.print = 5, pcs.compute=PCNUM, maxit = 500, weight.by.var = FALSE )
PCElbowPlot(object = EXP,num.pc=PCNUM)

PCAPlot(object = EXP, dim.1 = 1, dim.2 = 2)

PCUSE=1:10
EXP = RunTSNE(object = EXP, dims.use = PCUSE, do.fast = TRUE,check_duplicates = FALSE )

RES=0.6
EXP <- FindClusters(object = EXP, reduction.type = "pca", dims.use = PCUSE,  resolution = RES, print.output = 0, save.SNN = TRUE,force.recalc =T)

TSNEPlot(object = EXP,do.label=T)
