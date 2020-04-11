if (!require(ggplot2)) install.packages('ggplot2', repos='http://cran.us.r-project.org')
library(ggplot2)

if (!require(reshape2)) install.packages('reshape2', repos='http://cran.us.r-project.org')
library(reshape2)

if (!require(egg)) install.packages('egg', repos='http://cran.us.r-project.org')
library(egg)

args = commandArgs(trailingOnly=TRUE)

genes_coords_file = args[1]
coverage_per_window_file = args[2]
phased_snvs_file = args[3]
hap_blocks_file = args[4]
sv_coverage_file = args[5]
gene_coverage_file = args[6]
haplotype_alleles_file = args[7]

locus_phasing_plot = args[8]
sv_coverage_plot = args[9]
gene_coverage_plot = args[10]
haplotype_allele_plot = args[11]

genes = genes_coords_file
genes = read.table(genes)
genes$Region = substr(genes$V4,0,4)
genes$name = substr(genes$V4,4,100)
genes$pos = 1.5
genes$pos[c(1:nrow(genes)) %% 2 == 0] = -1

x = coverage_per_window_file
x = read.table(x)

y = phased_snvs_file
y = read.table(y)
y = subset(subset(y,V2 != "0/0"),V2 != "0/1")

z = hap_blocks_file
z = read.table(z)

p1 = ggplot() +
  geom_bar(data=x,aes(V2,V4,fill=as.character(V5)),stat="identity",position="stack") +
  theme_bw() + ylab("Read counts") + labs(fill="Haplotype") + xlab("IGH position") +
  scale_x_continuous(expand = c(0,0),limits=c(-1,1200000)) +
  scale_y_continuous(expand = c(0,0))

p2 = ggplot(data=genes,aes(x=V2,y=0,label=name,color=Region)) +
  geom_segment(data=genes,aes(x=V2,y=0,xend=V3,yend=0,size=10)) +
  theme_bw() + ylab("Genes") + scale_y_continuous(expand=c(0,0)) +
  geom_text(check_overlap = TRUE,size=3,angle=90,fontface = "bold") +
  scale_x_continuous(expand = c(0,0),limits=c(-1,1200000)) + guides(size=FALSE)

p3 = ggplot() +
  geom_point(data=y,aes(V1,V2,color=V2)) +
  theme_bw() + ylab("CCS SNVs") + labs(color="Genotype") +
  scale_x_continuous(expand = c(0,0),limits=c(-1,1200000))

p4 = ggplot() +
  geom_segment(data=z,aes(x=V1,y=V2,xend=V3,yend=1,size=5)) +
  theme_bw() + scale_y_continuous(expand = c(0,0)) + 
  ylab("Phased\nregions") + xlab("IGH position") + 
  scale_x_continuous(expand = c(0,0),limits=c(-1,1200000))

p1_2_3_4 = ggarrange(p1 + theme(axis.text.x = element_blank(), 
                                axis.ticks.x = element_blank(),
                                axis.title.x = element_blank()), 
                     p2 + theme(axis.text.x = element_blank(), 
                                axis.ticks.x = element_blank(),
                                axis.title.x = element_blank(),
                                axis.text.y = element_blank(),
                                axis.ticks.y = element_blank()),
                     p3 + theme(axis.text.x = element_blank(), 
                                axis.ticks.x = element_blank(),
                                axis.title.x = element_blank()),
                     p4 + theme(axis.text.y = element_blank(),
                                axis.ticks.y = element_blank(),
                                legend.position = "none"),
				nrow=4,heights=c(5,1,2,1))

pdf(locus_phasing_plot,width=14, height=10,units="in",res=150)
p1_2_3_4
dev.off()

sv_coverage = sv_coverage_file
sv_coverage = read.table(sv_coverage,header=TRUE)
sv_coverage = melt(sv_coverage, id.vars=c("chrom", "start","end","sv"))

sv_plot = ggplot(sv_coverage,aes(sv,value,fill=variable)) +
  geom_bar(stat="identity") +
  coord_flip() +
  theme_bw() + 
  xlab("Genotyped Structural\nVariant") + 
  ylab("Coverage") +
  labs(fill="Haplotype") +
  scale_fill_discrete(labels = c("0", "1", "2"))

png(sv_coverage_plot,width=6, height=2,units="in",res=300)
sv_plot
dev.off()

gene_coverage = gene_coverage_file
gene_coverage = read.table(gene_coverage,header=TRUE)
gene_coverage = melt(gene_coverage, id.vars=c("chrom", "start","end","gene"))
gene_coverage$region = substr(gene_coverage$gene,0,4)
gene_coverage$region_f = factor(gene_coverage$region, levels=c('IGHJ','IGHD','IGHV'))

gene_coverage = ggplot(gene_coverage,aes(reorder(gene,start),value,fill=variable)) +
  geom_bar(stat="identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,vjust=0.5,hjust=1)) +
  ylab("Coverage") + 
  xlab("Gene") +
  labs(fill="Haplotype") +
  scale_fill_discrete(labels = c("0", "1", "2")) +
  facet_grid(.~region_f,scales="free",space="free")

png(gene_coverage_plot,width=14, height=4,units="in",res=150)
gene_coverage
dev.off()

haplotype_alleles = haplotype_alleles_file
haplotype_alleles = read.table(haplotype_alleles,header=TRUE,sep="\t")
haplotype_alleles = melt(haplotype_alleles, id.vars=c("chrom", "start","end","gene_name"))
haplotype_alleles = subset(haplotype_alleles,variable == "haplotype_1_allele" | variable == "haplotype_2_allele")
haplotype_alleles_plot = ggplot(haplotype_alleles,aes(reorder(gene_name,start),variable)) +
  theme_bw() + 
  geom_tile(aes(fill = value),colour = "grey50") +
  theme(legend.position="bottom",axis.text.x = element_text(angle = 90, hjust = 1,vjust=.5),
        axis.text = element_text(colour = "black"),
        axis.title.x=element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y=element_blank(),
        strip.background = element_rect(linetype="solid",color="black", fill="gold3", size=1.5),
        panel.spacing = unit(0.1, "lines")) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  geom_point(data=subset(haplotype_alleles,value=="Deleted"),aes(gene_name,variable)) +
  geom_point(data=subset(haplotype_alleles,value=="Novel"),aes(gene_name,variable),shape=24) +
  guides(fill=guide_legend(ncol=20)) + labs(fill="Alleles")
png(haplotype_allele_plot,width=14, height=2,units="in",res=100)
haplotype_alleles_plot
dev.off()