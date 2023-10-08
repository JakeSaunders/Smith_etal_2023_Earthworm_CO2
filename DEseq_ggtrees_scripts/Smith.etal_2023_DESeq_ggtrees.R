# https://bioconductor.org/packages/release/bioc/html/ggtreeExtra.html
# https://bioconductor.org/packages/release/bioc/vignettes/ggtreeExtra/inst/doc/ggtreeExtra.html
# https://yulab-smu.top/treedata-book/chapter10.html

# Set dir and files names-----

dir <- "C:/Users/saundecj/Desktop/CO2 Earthworm/scripts/DEseq/"
setwd(dir)

xl.file <- "SmithEtal2023-canidates.xlsx"
ct.file <- "eh.1804.countTable.csv"

# read packages

# BiocManager::install("ggtreeExtra")
library(tidyverse)
library(DESeq2)
library(seqinr)
library(msa)
library(phyloseq)
library(ggthemes)
library(ggtree)
library(ggtreeExtra)
library(readxl)


# load data from trinotate and count table ----

# import trinotate candidates 
trinotate <- read_xlsx(xl.file) %>%
  mutate( transcript_id = sub(pattern = "TRINITY_DN","",transcript_id)) 

# import trinotate counttable from trinity 
counttable <- read.csv(ct.file)%>% 
  select(transcript_id, eh.head.S1.bam:eh.mid.S6.bam) %>% 
  mutate( transcript_id = sub(pattern = "TRINITY_DN","",transcript_id)) 

## Make Tidy candidates counttable for plotting  ----
counts <- inner_join( 
  x = counttable,
  y = trinotate %>% select(gene_family,transcript_id),
  by = "transcript_id") %>% 
  pivot_longer(
    cols = eh.head.S1.bam:eh.mid.S6.bam,
    names_to = c("eh","Tissue","Sample","File"),
    names_sep = "\\.",
    values_to = "count"
  ) %>% 
  select( transcript_id:gene_family,Tissue:Sample,count) %>% 
  group_by(transcript_id, gene_family, Tissue) 
  
counts$Tissue <- ifelse(counts$Tissue == "head","Prostomium","Midsegment")

counttable <- counttable %>% column_to_rownames(var = "transcript_id")

# DEseq Analysis  ----

# make df describing experimental design
coldata <- data.frame( 
  Tissue = factor(c(rep("Prostomium",3),rep("Midsegment",3)),levels = c("Prostomium","Midsegment")),
  row.names = colnames(counttable)
)

# check that sample names match
all(rownames(coldata) == colnames(counttable))

# load count table and experimental desing into DESeq 2 object
dds <- DESeqDataSetFromMatrix(countData = counttable,
                              colData = coldata,
                              design = ~ Tissue)

# run Differential expression analysis 
dds <- DESeq(dds)

# extract results from dds
res <- results(dds)

# reorder by p adjusted
res <- res[order(res$pvalue),]
summary(res)

# check if any of the significant genes are in out candidates
table( rownames(res[ which( res$padj < 0.05, res ),]) %in% trinotate$transcript_id)

# check significationly enriched for gene involved in CO2 detection
res.df <- as.data.frame(res)

res.df[ rownames(res) %in% trinotate$transcript_id, ]

# make data.frame of full trinotate data to work with
trinotate.full <- read_xlsx("eh.1804.trinotate_annotation_report.xlsx")
colnames(trinotate.full)[1] <- "gene_id"
trinotate.full <- trinotate.full %>% 
  mutate( transcript_id = sub(pattern = "TRINITY_DN","",transcript_id)) %>% 
  column_to_rownames(var = "transcript_id")

# subset 

rownames(res.df[ which( res.df$padj < 0.05, res.df ),] )



trinotate.full[ which( res.df$padj < 0.05, res.df ),]


trinotate.p05 <- trinotate.full[ rownames(trinotate.full) %in% rownames(res.df[ which( res.df$padj < 0.05, res.df ),] ), ]

write.csv(res.df[ which( res.df$padj < 0.05, res.df ),] ,file = "DEseq2.05.csv",row.names = T)
writexl::write_xlsx(trinotate.p05,path = "trinotate.p05.xlsx")

# remove large objects that are no longer needed 
rm(dds,res,coldata,trinotate.full)

# Make Tree figures ----

## Carbonic anhydrase ----
family <- "carbonic anhydrase"
title <-  "Carbonic Anhydrase"
cts <- counts %>% 
  filter(gene_family == family) %>% 
  ungroup() %>% 
  select(!gene_family) %>%
  group_by(transcript_id,Tissue) %>%
  summarise( Count =sum(count) ) %>% ungroup()

### calculate tree ----
seq <- AAStringSet(trinotate %>% filter( gene_family == family) %>% pull(peptide_seq))
names(seq) <- trinotate %>% filter( gene_family == family) %>% pull(transcript_id)
msa <- msaClustalW(inputSeqs = seq)
tree <- dist.alignment(msaConvert(msa, type="seqinr::alignment"),"identity")
tree <- ape::njs(tree)



### make df for annotations ----

anno <- left_join(
  x = trinotate %>% filter(gene_family == family) %>% 
    select(transcript_id, BLASTP.name,BLASTP.species), 
  y = read_csv("taxa_match.csv") %>% 
    select(-Taxa) %>% rename(BLASTP.species = Abbv),
  by = "BLASTP.species") %>% 
select(-BLASTP.species) %>% rename(Isoform = BLASTP.name)


anno$Isoform <- sub("Putative carbonic anhydrase-like protein","Carbonic anhydrase",anno$Isoform)
anno$Isoform <- sub(", mitochondrial|-related protein","",anno$Isoform)
anno[ anno$transcript_id == "76341_c1_g1_i9",]$Isoform <- "Carbonic anhydrase Z"
anno[ anno$transcript_id == "48383_c0_g1_i1",]$Isoform <- "Carbonic anhydrase 1"
anno[ is.na(anno$Isoform),]$Isoform <- "PF00194.20"
anno[ is.na(anno$Order),]$Order <- "Eukaryotic-type carbonic anhydrase"

anno$Isoform <- factor(anno$Isoform,levels=anno$Isoform)

# anno$Isoform <- factor(anno$Isoform, 
#   levels = c("Carbonic anhydrase 1","Carbonic anhydrase 2","Carbonic anhydrase 4",
#     "Carbonic anhydrase 5A","Carbonic anhydrase 6","Carbonic anhydrase 7","Carbonic anhydrase 9",
#     "Carbonic anhydrase 10","Carbonic anhydrase 12","Carbonic anhydrase 13","Carbonic anhydrase Z",
#     "Cell surface-binding protein","PF00194.20" ))

### draw tree ----
p <- ggtree(tree, layout="fan",branch.length = "none")
p <- rotate_tree(p, -90)

p <- p %<+% anno + 
  geom_tippoint(aes(shape=Isoform, color=Order),size=3,position = position_nudge(x=0.5))+
  scale_shape_manual(
    values = c(15:19,7:14,3), 
    guide=guide_legend(title = "Gene Annotation", keywidth = 0.5, keyheight = 0.5, 
      order=1,override.aes=list(shape=c(15:19,7:14,3)))
  ) + 
  scale_color_manual(
    values = c("#ffa600","#2f4b7c","Black","#a05195","#d45087","#f95d6a","#ff7c43","#003f5c"), 
    guide=guide_legend(title = "Order of BLASTP Match", keywidth = 0.5, keyheight = 0.5, order=2)
  ) 

p +geom_fruit(
    data= cts,
    geom=geom_col,
    mapping=aes( x= Count , y=transcript_id, fill = Tissue),
    orientation="y",
    position = position_dodgex(),
    offset = 0.2,
    axis.params=list(
       axis       = "x",
       text.size  = 1.8,
       nbreak     = 3,
     ),
  ) +
  scale_fill_manual(values=c("dodgerblue","darkgreen"),
    guide=guide_legend(title = "Tissue Transcript Counts", 
    keywidth = 0.5, keyheight = 0.5, order=2)
  ) +
  geom_tiplab(size=3, color="grey40",hjust = -0.75) +
  labs(title = paste( title, "Transcripts") ) +
  hexpand(.2, direction = -1)

ggsave(paste0("Tree - ",title, ".pdf"),height = 8,width = 8,units = "in")
### clean up environment ----
rm(anno,cts,msa,p,seq,tree,family,title)
## ASIC ----
family <- "ASIC"
title <-  "Acid Sensing Ion Channels"
cts <- counts %>% 
  filter(gene_family == family) %>% 
  ungroup() %>% 
  select(!gene_family) %>%
  group_by(transcript_id,Tissue) %>%
  summarise( Count =sum(count) ) %>% ungroup()

### calculate tree ----
seq <- AAStringSet(trinotate %>% filter( gene_family == family) %>% pull(peptide_seq))
names(seq) <- trinotate %>% filter( gene_family == family) %>% pull(transcript_id)
msa <- msaClustalW(inputSeqs = seq)
tree <- dist.alignment(msaConvert(msa, type="seqinr::alignment"),"identity")
tree <- ape::njs(tree)

### make df for annotations ----

anno <- left_join(
  x = trinotate %>% filter(gene_family == family) %>% 
    select(transcript_id, BLASTP.name,BLASTP.species), 
  y = read_csv("taxa_match.csv") %>% 
    select(-Taxa) %>% rename(BLASTP.species = Abbv),
  by = "BLASTP.species") %>% 
  select(-BLASTP.species) %>% rename(Isoform = BLASTP.name)


anno$Isoform <- sub(" \\{ECO\\:0000303\\|PubMed\\:9062189\\}","",anno$Isoform)
anno$Isoform <- sub(" \\{ECO\\:0000303\\|PubMed\\:10798398\\}","",anno$Isoform)
anno$Isoform <- sub(" \\{ECO\\:0000303\\|PubMed\\:10842183\\}","",anno$Isoform)

#anno$Isoform <- factor(anno$Isoform,levels=sort(unique(anno$Isoform)))

### draw tree ----
p <- ggtree(tree, layout="fan",branch.length = "none")
p <- rotate_tree(p, -90)

p <- p %<+% anno + 
  geom_tippoint(aes(shape=Isoform, color=Order),size=3,position = position_nudge(x=0.5))+
  scale_shape_manual(
    values = c(15:18,8:10,13), 
    guide=guide_legend(title = "Gene Annotation", keywidth = 0.5, keyheight = 0.5, 
                       order=1,override.aes=list(shape=c(15:18,8:10,13)))
  ) +
  scale_color_manual(
    values = c("#003f5c","#58508d","#bc5090","#ff6361","#ffa600"), 
    guide=guide_legend(title = "Order of BLASTP Match", keywidth = 0.5, keyheight = 0.5, order=2)
  ) 


p +geom_fruit(
  data= cts,
  geom=geom_col,
  mapping=aes( x= Count , y=transcript_id, fill = Tissue),
  orientation="y",
  position = position_dodgex(),
  offset = 0.2,
  axis.params=list(
    axis       = "x",
    text.size  = 1.8,
    nbreak     = 3,
  ),
) +
  scale_fill_manual(values=c("dodgerblue","darkgreen"),
                    guide=guide_legend(title = "Tissue Transcript Counts", 
                                       keywidth = 0.5, keyheight = 0.5, order=2)
  ) +
  geom_tiplab(size=3, color="grey40",hjust = -0.75) +
  labs(title = paste( title, "Transcripts") ) +
  hexpand(.2, direction = -1)

ggsave(paste0("Tree - ",title, ".pdf"),height = 8,width = 8,units = "in")

### clean up environment ----
rm(anno,cts,msa,p,seq,tree,family,title)

## guanylate cyclase ----
family <- "guanylate cyclase"
title <-  "Guanylate Cyclase"
cts <- counts %>% 
  filter(gene_family == family) %>% 
  ungroup() %>% 
  select(!gene_family) %>%
  group_by(transcript_id,Tissue) %>%
  summarise( Count =sum(count) ) %>% ungroup()

### calculate tree ----
seq <- AAStringSet(trinotate %>% filter( gene_family == family) %>% pull(peptide_seq))
names(seq) <- trinotate %>% filter( gene_family == family) %>% pull(transcript_id)
msa <- msaClustalW(inputSeqs = seq)
tree <- dist.alignment(msaConvert(msa, type="seqinr::alignment"),"identity")
tree <- ape::njs(tree)

### make df for annotations ----

anno <- left_join(
  x = trinotate %>% filter(gene_family == family) %>% 
    select(transcript_id, BLASTP.name,BLASTP.species), 
  y = read_csv("taxa_match.csv") %>% 
    select(-Taxa) %>% rename(BLASTP.species = Abbv),
  by = "BLASTP.species") %>% 
  select(-BLASTP.species) %>% rename(Isoform = BLASTP.name)

anno$Isoform <- sub(" \\{ECO\\:0000305\\}","",anno$Isoform)

unique(anno$Isoform)
unique(anno$Order)

#anno$Isoform <- factor(anno$Isoform,levels=sort(unique(anno$Isoform)))

### draw tree ----
p <- ggtree(tree, layout="fan",branch.length = "none")
p <- rotate_tree(p, -90)

p <- p %<+% anno + 
  geom_tippoint(aes(shape=Isoform, color=Order),size=3,position = position_nudge(x=0.5))+
  scale_shape_manual(
    values = c(15:18,8:10), 
    guide=guide_legend(title = "Gene Annotation", keywidth = 0.5, keyheight = 0.5, 
                       order=1,override.aes=list(shape=c(15:18,8:10)))
  ) +
  scale_color_manual(
    values = c("#ffa600","#7a5195","#ef5675","#003f5c"), 
    guide=guide_legend(title = "Order of BLASTP Match", keywidth = 0.5, keyheight = 0.5, order=2)
  ) 

p +geom_fruit(
  data= cts,
  geom=geom_col,
  mapping=aes( x= Count , y=transcript_id, fill = Tissue),
  orientation="y",
  position = position_dodgex(),
  offset = 0.2,
  axis.params=list(
    axis       = "x",
    text.size  = 1.8,
    nbreak     = 3,
  ),
) +
  scale_fill_manual(values=c("dodgerblue","darkgreen"),
                    guide=guide_legend(title = "Tissue Transcript Counts", 
                                       keywidth = 0.5, keyheight = 0.5, order=3)
  ) +
  geom_tiplab(size=3, color="grey40",hjust = -0.75) +
  labs(title = paste( title, "Transcripts") ) +
  hexpand(.2, direction = -1)

ggsave(paste0("Tree - ",title, ".pdf"),height = 8,width = 8,units = "in")

### clean up environment ----
rm(anno,cts,msa,p,seq,tree,family,title)

## Otopetrin ----
family <- "Otopetrin"
title <-  "Otopetrin"
cts <- counts %>% 
  filter(gene_family == family) %>% 
  ungroup() %>% 
  select(!gene_family) %>%
  group_by(transcript_id,Tissue) %>%
  summarise( Count =sum(count) ) %>% ungroup()

### calculate tree ----
seq <- AAStringSet(trinotate %>% filter( gene_family == family) %>% pull(peptide_seq))
names(seq) <- trinotate %>% filter( gene_family == family) %>% pull(transcript_id)
msa <- msaClustalW(inputSeqs = seq)
tree <- dist.alignment(msaConvert(msa, type="seqinr::alignment"),"identity")
tree <- ape::njs(tree)

### make df for annotations ----

anno <- left_join(
  x = trinotate %>% filter(gene_family == family) %>% 
    select(transcript_id, BLASTP.name,BLASTP.species), 
  y = read_csv("taxa_match.csv") %>% 
    select(-Taxa) %>% rename(BLASTP.species = Abbv),
  by = "BLASTP.species") %>% 
  select(-BLASTP.species) %>% rename(Isoform = BLASTP.name)


anno$Isoform <- ifelse(is.na(anno$Isoform),yes = "PF03189.12",no="OtopLc")
anno$Order <- ifelse(is.na(anno$Order),yes = "Otopetrin Protein Domain",no=anno$Order)

unique(anno$Isoform)
unique(anno$Order)

#anno$Isoform <- factor(anno$Isoform,levels=sort(unique(anno$Isoform)))

### draw tree ----
p <- ggtree(tree, layout="fan",branch.length = "none")
p <- rotate_tree(p, -90)

p %<+% anno + 
  geom_tippoint(aes(shape=Isoform, color=Order),size=3,position = position_nudge(x=0.5))

p <- p %<+% anno + 
  geom_tippoint(aes(shape=Isoform, color=Order),size=3,position = position_nudge(x=0.5))+
  scale_shape_manual(
    values = c(16,17), 
    guide=guide_legend(title = "Gene Annotation", keywidth = 0.5, keyheight = 0.5, 
                       order=1,override.aes=list(shape=c(16,17)))
  ) +
  scale_color_manual(
    values = c("#ffa600","black","#ef5675","#003f5c"), 
    guide=guide_legend(title = "Order of BLASTP Match", keywidth = 0.5, keyheight = 0.5, order=2)
  ) 

p +geom_fruit(
  data= cts,
  geom=geom_col,
  mapping=aes( x= Count , y=transcript_id, fill = Tissue),
  orientation="y",
  position = position_dodgex(),
  offset = 0.2,
  axis.params=list(
    axis       = "x",
    text.size  = 1.8,
    nbreak     = 3,
  ),
) +
  scale_fill_manual(values=c("dodgerblue","darkgreen"),
                    guide=guide_legend(title = "Tissue Transcript Counts", 
                                       keywidth = 0.5, keyheight = 0.5, order=3)
  ) +
  geom_tiplab(size=3, color="grey40",hjust = -0.75) +
  labs(title = paste( title, "Transcripts") ) +
  hexpand(.2, direction = -1)

ggsave(paste0("Tree - ",title, ".pdf"),height = 8,width = 8,units = "in")

### clean up environment ----
rm(anno,cts,msa,p,seq,tree,family,title)
## TRPA ----

family <- "TRPA"
title <-  "TRPA"
cts <- counts %>% 
  filter(gene_family == family) %>% 
  ungroup() %>% 
  select(!gene_family) %>%
  group_by(transcript_id,Tissue) %>%
  summarise( Count =sum(count) ) %>% ungroup()

### calculate tree ----
seq <- AAStringSet(trinotate %>% filter( gene_family == family) %>% pull(peptide_seq))
names(seq) <- trinotate %>% filter( gene_family == family) %>% pull(transcript_id)
msa <- msaClustalW(inputSeqs = seq)
tree <- dist.alignment(msaConvert(msa, type="seqinr::alignment"),"identity")
tree <- ape::njs(tree)

### make df for annotations ----

anno <- left_join(
  x = trinotate %>% filter(gene_family == family) %>% 
    select(transcript_id, BLASTP.name,BLASTP.species), 
  y = read_csv("taxa_match.csv") %>% 
    select(-Taxa) %>% rename(BLASTP.species = Abbv),
  by = "BLASTP.species") %>% 
  select(-BLASTP.species) %>% rename(Isoform = BLASTP.name)


anno$Isoform <- sub(
  pattern = "Transient receptor potential cation channel subfamily A member 1",
  replacement = "TRPA1",x = anno$Isoform)

unique(anno$Isoform)
unique(anno$Order)

#anno$Isoform <- factor(anno$Isoform,levels=sort(unique(anno$Isoform)))

### draw tree ----
p <- ggtree(tree, layout="fan",branch.length = "none")
p <- rotate_tree(p, -90)

p %<+% anno + 
  geom_tippoint(aes(shape=Isoform, color=Order),size=3,position = position_nudge(x=0.5))

p <- p %<+% anno + 
  geom_tippoint(aes(shape=Isoform, color=Order),size=3,position = position_nudge(x=0.2))+
  scale_shape_manual(
    values = c(16), 
    guide=guide_legend(title = "Gene Annotation", keywidth = 0.5, keyheight = 0.5, 
                       order=1,override.aes=list(shape=c(16)))
  ) +
  scale_color_manual(
    values = c("#ffa600","#bc5090","#003f5c"), 
    guide=guide_legend(title = "Order of BLASTP Match", keywidth = 0.5, keyheight = 0.5, order=2)
  ) 

p +geom_fruit(
  data= cts,
  geom=geom_col,
  mapping=aes( x= Count , y=transcript_id, fill = Tissue),
  orientation="y",
  position = position_dodgex(),
  offset = 0.2,
  axis.params=list(
    axis       = "x",
    text.size  = 1.8,
    nbreak     = 3,
  ),
) +
  scale_fill_manual(values=c("dodgerblue","darkgreen"),
                    guide=guide_legend(title = "Tissue Transcript Counts", 
                                       keywidth = 0.5, keyheight = 0.5, order=3)
  ) +
  geom_tiplab(size=3, color="grey40",hjust = -0.75) +
  labs(title = paste( title, "Transcripts") ) +
  hexpand(.2, direction = -1)

ggsave(paste0("Tree - ",title, ".pdf"),height = 8,width = 8,units = "in")

### clean up environment ----
rm(anno,cts,msa,p,seq,tree,family,title)
  