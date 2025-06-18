library(microbiome)
library(knitr);library(phyloseq);require(dplyr)
psd=readRDS("C:/TESIS/ps_Data_ash")
sdata <- read.csv("C:/TESIS/table data 16S.csv",row.names = 1)

psd <- prune_samples(row.names(sdata),x = psd)
reads_sample <- readcount(psd)
x=as.data.frame(sample_data(psd))
sample_data(psd) <- sdata
#filtrado de prevalencia de taxones a partir de al menos el 20% del total de las muestras
psd <- filter_taxa(psd, function(x) sum(x > 2) > (0.2*length(x)), TRUE)
psd1 = prune_samples(sample_sums(psd) > 1000, psd)
prevdf = apply(X = otu_table(psd1),
               MARGIN = ifelse(taxa_are_rows(psd1), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Le agregamos la taxonomía
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(psd1),
                    tax_table(psd))
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))}) -> dfprev
kable(dfprev)
#remover taxa que no corresponde a microorganismos como cloroplastos, mitocondrias y otros
filterPhyla2 <- c("Chloroplast", "Mitochondria", "Eukaryota")
psd1 <- subset_taxa(psd1, !Kingdom %in% filterPhyla2)
psd1 <- subset_taxa(psd1, !Phylum %in% filterPhyla2)
psd1 <- subset_taxa(psd1, !Class %in% filterPhyla2)
psd1 <- subset_taxa(psd1, !Order %in% filterPhyla2)
psd1 <- subset_taxa(psd1, !Family %in% filterPhyla2)
psd1 <- subset_taxa(psd1, !Genus %in% filterPhyla2)
#seleccionar taxas de interes
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(psd1, "Phylum"))
library(ggplot2)
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(psd1),color=Phylum)) +
# Agregamos una línea para nuestro umbral
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")
# Definimos el umbral de prevalencia a un 5%
(prevalenceThreshold = 0.05 * nsamples(psd1))
# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
(psd2 = prune_taxa(keepTaxa, psd1))
reads_sample2 <- readcount(psd2)
# Reemplazamos las secuencias por un nombre genérico
taxa_names(psd2) <- paste0("ASV", seq(ntaxa(psd2)))
sample_sum_df <- data.frame(sum = sample_sums(psd2))
### sacar tabla 1

table_QC <- data.frame(BEFORE_FILTER =reads_sample,AFTER_FILTER=reads_sample2,Q30=NA)# Q30 se pega de lo de NOVOGENE
write.csv(table_QC,"C:/TESIS/QC.csv",na = "",row.names = T)
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "purple", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank()) 
# Primero cargamos algunos scripts de manera remota
scripts <- c("graphical_methods.R",
             "tree_methods.R",
             "plot_merged_trees.R",
             "specificity_methods.R",
             "ternary_plot.R",
             "richness.R",
             "edgePCA.R",
             "copy_number_correction.R",
             "import_frogs.R",
             "prevalence.R",
             "compute_niche.R")
urls <- paste0("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/R/", scripts)

for (url in urls) {
  source(url)
}
#p <- ggrare(psd2, step = 100, color = "species", label = "sample_ID", se = TRUE)
#(p <- p + facet_wrap(~species))
#p <- ggrare(psd2, step = 100, color = "phylum", label = "sample_ID", se = TRUE)
#(p <- p + facet_wrap(~phylum))
#install.packages("BiocManager")
#BiocManager::install("microbiome")
library(microbiome)
#sample_data(psd2)$group=factor(sample_data(psd2)$group)
#res <- plot_frequencies(sample_data(psd2), "group", "phylum")
#print(res$plot)
# Transformamos conteos en porcentaje
psd2r  = transform_sample_counts(psd2, function(x) x / sum(x) )

# Filtramos las taxa con una abundancia inferior al 1%
#(psd2r.filtrado = filter_taxa(psd2r, function(x) sum(x) > 1, TRUE))
plot_richness(psd2, x = "group", 
              measures = c("Observed", "Chao1", "Shannon","Simpson"))  + 
  geom_boxplot(alpha=.7) 
  #scale_color_manual(values = c("#a6cee3", "#b2df8a", "#fdbf6f")) + 
  #scale_fill_manual(values = c("#a6cee3", "#b2df8a", "#fdbf6f"))
plot_richness(psd2, x = "line", 
              measures = c("Observed", "Chao1", "Shannon","Simpson"))  + 
  geom_boxplot(alpha=.7) 
# Guardamos un dataframe con las medidas de diversidad alfa

#devtools::install_github('twbattaglia/btools')
library(btools)
#alpha_pd <- estimate_pd(psd2)
# Combinamos la metadata con alpha.diversity
#data <- cbind(sample_data(psd2), alpha_pd) 
####################################################################

otu_jac <- as.data.frame(psd2@otu_table)
otu_jac[otu_jac>0] <- 1
# Y calculamos un ANOVA
#psd2.anova <- aov(PD ~ species, data) 
# install.packages("xtable")
library(xtable)
tab <- global(psd2, index = "all")
require(magrittr)
#write.csv()
#Shapiro-Wilk normality test
#par(mfrow = c(1, 3))

#Then plot each metric
hist(tab$diversity_shannon, main="Shannon", xlab="", breaks=10)
hist(tab$dominance_simpson, main="Simpson", xlab="", breaks=10)
hist(tab$chao1, main="Chao1", xlab="", breaks=15)
#hist(tab$, main="ACE richness", xlab="", breaks=15)


require(kableExtra)
head(tab) %>%
  kable(format = "html", col.names = colnames(tab), digits = 2) %>%
  kable_styling() %>%
  kableExtra::scroll_box(width = "100%", height = "300px")

#data <- cbind(sample_data(psd2), tab) 
#psd2.anova <- aov(diversity_shannon~ group, data) 

#summary(psd2.anova)

x_k <- kruskal.test(tab$observed ~ sdata$year)
x_k$p.value


x_k <- kruskal.test(tab$observed ~ sdata$line)
x_k$p.value

x_k <- kruskal.test(tab$diversity_shannon ~ sdata$line)
x_k$p.value

x_k <- kruskal.test(tab$diversity_gini_simpson ~ sdata$line)
x_k$p.value

x_k <- wilcox.test(tab$observed ~ sdata$line)
x_k$p.value

x_k <- kruskal.test(tab$observed ~ sdata$pedigree)
x_k$p.value

x_k <- kruskal.test(tab$chao1 ~ sdata$group)
x_k$p.value

x_k <- wilcox.test(tab$chao1 ~ sdata$line)
x_k$p.value

x_k <- kruskal.test(tab$chao1 ~ sdata$pedigree)
x_k$p.value

x_k <- kruskal.test(tab$diversity_shannon ~ sdata$group)
x_k$p.value

x_k <- kruskal.test(tab$diversity_shannon ~ data$pedigree)
x_k$p.value

x_k <- kruskal.test(tab$diversity_shannon ~ sdata$year)
x_k$p.value
#si son dos grupos es wilcox sino es krusktal
x_k <- wilcox.test(tab$dominance_simpson ~ sdata$line)
x_k$p.value

x_k <- wilcox.test(tab$dominance_simpson ~ sdata$pedigree)
x_k$p.value

p4cols <- c("blue","lightblue","red","orange")

psd2.dpcoa_ori<- ordinate(psd2r, method = "PCoA", distance = "bray")
evals <- psd2.dpcoa_ori$eig
pord3 <- plot_ordination(psd2r, psd2.dpcoa_ori, color = "group") +
  labs(col = "group") +
  geom_point(size=5)

pord3 = pord3 + scale_color_manual(values = p4cols)

pord3 + facet_grid(. ~ line)#,cols = 2,rows=2)
pord3
#psd2.dpcoa <- ordinate(psd2r, method = "NMDS", distance = "bray")
#evals <- psd2.dpcoa$eig
#pord3 <- plot_ordination(year, psd2.dpcoa, color = "group", shape = "group") +
  #labs(col = "group") +
  #coord_fixed(sqrt(evals[2] / evals[1])) +
  #scale_color_manual(values = c("#a6cee3", "#b2df8a", "#fdbf6f")) + 
  #scale_fill_manual(values = c("#a6cee3", "#b2df8a", "#fdbf6f")) +
  #geom_point(size=4)
#pord3
year=subset_samples(psd2r, year=="2021")
psd2.dpcoa <- ordinate(year, method = "PCoA", distance = "bray")
evals <- psd2.dpcoa$eig
pord3 <- plot_ordination(year, psd2.dpcoa, color = "group", shape = "group") +
  labs(col = "group") +
    geom_point(size=5)
pord3
year=subset_samples(psd2r, year=="2021")
psd2.dpcoa <- ordinate(year, method = "PCoA", distance = "bray")
evals <- psd2.dpcoa$eig
pord3 <- plot_ordination(year, psd2.dpcoa, color = "group",
                        shape = "group") +
  labs(col = "group") +
  geom_point(size=4)
#pord3

gg_color_hue <- function(n){
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}
#p4cols <- gg_color_hue(length(unique(year@sam_data$group)))
p4cols <- c("blue","red")
pord3 = pord3 + scale_color_manual(values = p4cols)
pord3
## agglomerate at the Family taxonomic rank
(x1 <- tax_glom(psd2, taxrank="Family") )
#BiocManager::install(c("phyloseq", "microbiome", "ComplexHeatmap"), update = FALSE)
##install.packages(
  #"microViz",
  #repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
#)
#install.packages("ggtext") # for rotated labels on ord_plot() 
#install.packages("ggraph") # for taxatree_plots()
#install.packages("DT") # for tax_fix_interactive()
#install.packages("corncob") # for example datasets and beta binomial models
# illustrative simple customised example
 require(microViz)
comp_barplot(psd2,
    tax_level = "Family", n_taxa = 8,
    bar_outline_colour = NA,
    sample_order = "bray",
    bar_width = 0.7,
    taxon_renamer = toupper
  ) + coord_flip()
#> Registered S3 method overwritten by 'seriation':
#>   method         from 
#>   reorder.hclust vegan
# Necesitamos obtener las taxa más abundantes, e.g top 15
#devtools::install_github("gmteunisse/fantaxtic")
# Ya que no todas las taxa fueron clasificadas a nivel de especie, generamos etiquetas compuestas de distintos rangos taxonómicos para el gráfico
#top15 <- nested_top_taxa(psd2, n_top_taxa = 15,top_tax_level= "Kingdom", nested_tax_level= "Species",by_proportion= TRUE)
#top15 <- fantaxtic::top_taxa(psd2, n_taxa = 15)#, #relative = T,
#library(fantaxtic)
#top15 <- name_na_taxa(top15$ps_obj, 
                      #include_rank = T)
#plot_nested_bar(ps_obj = top15,
                #top_level = "Phylum",
               # nested_level = "Species",
                #include_rank = T) +
#facet_wrap(~group,
          # scales = "free_x") +
  #labs(title = "") +
  #theme(plot.title = element_text(hjust = 0.5, 
                                 # size = 8, 
                                  #face = "bold"),
        #legend.key.size = unit(10, 
                               #"points"))


require(microViz)
psd2_copy <- psd2
psd2_copy <- psd2_copy %>% tax_fix()
# make a list of 2 harmonised composition plots (grouped by sex)
p <- comp_barplot(psd2_copy,
                  n_taxa = 8, tax_level = "Phylum",
                  bar_outline_colour = "black", merge_other = TRUE,
                  #sample_order = "aitchison",
                  group_by = "group"
)

# plot them side by side with patchwork package
patch <- patchwork::wrap_plots(p, ncol = 2, guides = "collect")
patch & coord_flip() # make bars in all plots horizontal (note: use & instead of +)
#discard_other = T, other_label = "Other")
# Ya que no todas las taxa fueron clasificadas a nivel de especie, generamos etiquetas compuestas de distintos rangos taxonómicos para el gráfico

# Finalmente graficamos

#require(fantaxtic)
ptop15= fantaxtic::fantaxtic_bar(top15, 
                                 color_by = "Phylum",
                                 label_by = "Phylum",
                                 facet_by = "line", 
                                 #grid_by = "pedigree",
                                 order_alg ="hclust",
                                 other_color = "Grey",
                                 facet_cols = 2,
                                 gen_uniq_lbls = TRUE)
#ptop15
#install.packages("ampvis2")
#remotes::install_github("kasperskytte/ampvis2")
library(ampvis2)
obj <- psd2
# Cambiamos la orientación de la otu_table
t(otu_table(obj)) -> otu_table(obj)
# Extraemos las tablas
otutable <- data.frame(OTU = rownames(phyloseq::otu_table(obj)@.Data),
                       phyloseq::otu_table(obj)@.Data,
                       phyloseq::tax_table(obj)@.Data,
                       check.names = FALSE)
# Extraemos la metadada
metadata <- data.frame(phyloseq::sample_data(obj), 
                       check.names = FALSE)
metadata$SampleID = rownames(metadata)
metadata= metadata[,c(7,1:6)]
# ampvis2 requiere que 1) los rangos taxonómicos sean siete y vayan de Kingdom a Species y 2) la primera columna de la metadata sea el identificador de cada muestra
# Entonces duplicamos la columna Género y le cambiamos el nombre a Species
#otutable$Species = otutable$Genus
# Reordenamos la metadata
#metadata <- metadata[,c("run","sample_ID","bioproject_accession","study","biosample_accession","experiment","SRA_Sample","geo_loc_name","collection_date","sample_type","species","common_name","AvgSpotLen")]

# finalmente generamos el objeto ampvis
av2 <- amp_load(otutable, metadata)
amp_heatmap(av2, 
            group_by = "year", 
            facet_by = "genotype", 
            plot_values = TRUE,
            tax_show = 20,
            tax_aggregate = "Genus",
            tax_add = "Phylum",
            plot_colorscale = "sqrt",
            plot_legendbreaks = c(1, 5, 10)
)
amp_heatmap(av2, 
            group_by = "year", 
            facet_by = "group", 
            plot_values = TRUE,
            tax_show = 20,
            tax_aggregate = "Genus",
            tax_add = "Phylum",
            plot_colorscale = "sqrt",
            plot_legendbreaks = c(1, 5, 10))
amp_heatmap(av2, 
            group_by = "group", 
            facet_by = "line", 
            plot_values = TRUE,
            tax_show = 20,
            tax_aggregate = "Phylum",
            tax_add = "Phylum",
            plot_colorscale = "sqrt",
            plot_legendbreaks = c(1, 5, 10))
amp_heatmap(av2, 
            group_by = "group", 
            facet_by = "line", 
            plot_values = TRUE,
            tax_show = 20,
            tax_aggregate = "Genus",
            tax_add = "Phylum",
            plot_colorscale = "sqrt",
            plot_legendbreaks = c(1, 5, 10))
#amp_venn(av2, group_by = "", cut_a = 0, cut_f = 50, text_size = 3)
#Kruskal-Wallis rank sum test
require(agricolae)
#install.packages("agricolae")
otu_table <- psd2r@otu_table

require(parallel)
krus_Df <- data.frame(matrix(ncol = 3,nrow = ncol(otu_table)))
colnames(krus_Df)[1:3] <- c("ASV","KRUSKAL","FDR")
parallel::detectCores()
cl <- parallel::makeCluster(6)  #si son 8 cores nunca se usan 8, maximo 7!
parallel::clusterExport(cl, varlist=c("metadata","otu_table"),envir=environment())
xx_krus <- parallel::parLapplyLB(cl,
                                   X = seq_len(ncol(otu_table)),
                                   fun = function (i){

              x_k <- kruskal.test(otu_table[,i] ~ metadata$line)
              x_k$p.value
})
parallel::stopCluster(cl)
krus_Df$ASV <- colnames(otu_table)
krus_Df$KRUSKAL <-  unlist(xx_krus)
krus_Df$FDR <- p.adjust(krus_Df$KRUSKAL,method = "fdr")
#organizar de menor a mayor
krus_Df <- krus_Df[order(krus_Df$FDR,decreasing = F),]
# la unica cercana a 0.05 ()
#psd2r@tax_table[3177,]
####
tax_table <- as.data.frame(psd2r@tax_table)
tax_table$ASV <- row.names(tax_table)
tax_table <- tax_table[,c(8,1:7)]
krus_Df_filtered <- left_join(krus_Df,tax_table,"ASV")

################################################################################
#group
####
krus_Df2 <- data.frame(matrix(ncol = 3,nrow = ncol(otu_table)))
colnames(krus_Df2)[1:3] <- c("ASV","KRUSKAL","FDR")
parallel::detectCores()
cl <- parallel::makeCluster(6)  #si son 8 cores nunca se usan 8, maximo 7!
parallel::clusterExport(cl, varlist=c("metadata","otu_table"),envir=environment())
xx_krus2 <- parallel::parLapplyLB(cl,
                                 X = seq_len(ncol(otu_table)),
                                 fun = function (i){
                                   
                                   x_k <- kruskal.test(otu_table[,i] ~ metadata$group)
                                   x_k$p.value
                                 })
parallel::stopCluster(cl)
krus_Df2$ASV <- colnames(otu_table)
krus_Df2$KRUSKAL <-  unlist(xx_krus2)
krus_Df2$FDR <- p.adjust(krus_Df2$KRUSKAL,method = "fdr")
#organizar de menor a mayor
krus_Df2 <- krus_Df2[order(krus_Df2$FDR,decreasing = F),]

krus_Df2_filtered <- krus_Df2[which(krus_Df2$FDR<0.05),]
####
tax_table <- as.data.frame(psd2r@tax_table)
tax_table$ASV <- row.names(tax_table)
tax_table <- tax_table[,c(8,1:7)]
krus_Df2_filtered <- left_join(krus_Df2_filtered,tax_table,"ASV")

##conteo
tapply(krus_Df2_filtered$Phylum,krus_Df2_filtered$Phylum,length)
xk_list <- list()
for(i in 1:nrow(krus_Df2_filtered)){
  xk <- kruskal(otu_table[,krus_Df2_filtered$ASV[i]],metadata$group,group=F,p.adj="fdr")
  xk <- xk$comparison
  xk <- xk[which(xk$Signif.!=" "),]
  xk$ASV <- krus_Df2_filtered$ASV[i]
  xk$groups <- row.names(xk)
  xk <- xk[,c(4,5,1,2,3)]
  xk_list[[i]] <- xk
}
xk_list <- do.call(rbind,xk_list)
xk_list2 <- left_join(xk_list,tax_table,"ASV")

#guardar kruskal
write.csv(krus_Df2_filtered,"C:/TESIS/kruskal_groups2.csv",row.names = F, quote = F,na = "")
write.csv(xk_list2,"C:/TESIS/kruskal_groups_wilcoxon.csv",row.names = F, quote = F,na = "")

####guardar tax, otu, etc..
write.csv(tax_table,"C:/TESIS/tax_table.csv",row.names = F, quote = F,na = "")
write.csv(otu_table,"C:/TESIS/otu_table_trans.csv",row.names = T, quote = F,na = "")
write.csv(psd2@otu_table,"C:/TESIS/otu_table.csv",row.names = T, quote = F,na = "")

#Calculate distance and save as a matrix
BC.dist=vegdist(otu_table(psd2r), distance="bray")
#Run PERMANOVA on distances
disp.age = betadisper(BC.dist, metadata$group)
permutest(disp.age, pairwise=F, permutations=10000)
#perm <- anosim(BC.dist, metadata$group, permutations = 1000)

perm <- adonis2(BC.dist ~ metadata$group, permutations = 10000)
perm
#line
disp.age2 = betadisper(BC.dist, metadata$line)
permutest(disp.age2, pairwise=F, permutations=1000)
perm2 <- adonis2(BC.dist ~ metadata$line, permutations = 10000)
perm2
#year
disp.age2 = betadisper(BC.dist, metadata$year)
permutest(disp.age2, pairwise=F, permutations=10000)
perm2 <- adonis2(BC.dist ~ metadata$year, permutations = 10000)
perm2
#line and heat interaction
perm2 <- adonis2(BC.dist ~ metadata$line*metadata$group, permutations = 10000)
perm2
#line and year interaction
perm2 <- adonis2(BC.dist ~ metadata$line*metadata$year, permutations = 10000)
perm2

#heat and wheat interaction
perm2 <- adonis2(BC.dist ~ metadata$group*metadata$line, permutations = 10000)
perm2

# 
# x <- otu_table(psd2r)
# x <- x[c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48)]
# m <- metadata$line[c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48)]
# BC.dist2=vegdist(x, distance="bray")
# perm2 <- adonis2(BC.dist2 ~ m,
#                  permutations = 10000)
# perm2
# x <- otu_table(psd2r)
# x <- x[c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48)-1]
# m <- metadata$line[c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48)-1]
# perm2 <- adonis2(BC.dist2 ~ m,
#                  permutations = 10000)
# perm2
################################################################################
cor.spearman = as.data.frame(cor(psd2r@otu_table, metadata$cuant.canoppy.temperature, method = "spearman"))
cor.spearman$ASV <- row.names(cor.spearman)
cor.spearman2 <- cor.spearman[which(cor.spearman$V1 > 0.7),]
cor.spearman3 <- cor.spearman[which(cor.spearman$V1 < -0.7),]

cor.spearman <- rbind(cor.spearman2,cor.spearman3) 
cor.spearman <-cor.spearman[,c(2,1)]
#cor.spearman <- 
tt <- as.data.frame(tax_table(psd2r))
tt$ASV <- row.names(tt)
cor.spearman <- left_join(cor.spearman,tt,"ASV")
write.csv(cor.spearman,"C:/TESIS/correlation_cannopy_rel.csv",row.names = F, quote = F,na = "")
###


library(UpSetR)
library(ggplot2)
#ggVennDiagram::ggVennDiagram(x) + scale_fill_gradient(low="grey90",high = "red")

x <- as.data.frame(psd2@otu_table)
#https://www.yanh.org/2021/01/01/microbiome-r/
#glom <- tax_glom(psd2, taxrank = 'Family', NArm = T)
#x <- as.data.frame(glom@otu_table)
#colnames(x) <- glom@tax_table[,5]#5 family, 6 genus
#x$line <- metadata$line
#x$line <- paste0(metadata$line,"-",metadata$group)
#x$line <- paste0(glom@sam_data$line,"-",glom@sam_data$group)
x$line <- metadata$group
lines <- unique(x$line)
#line
#line group
#lines <- lines[c(1,2,3,4,5,6,7,8)]
#H 2021
#lines <- lines[c(1,3,5,7)]
#H 2022
#lines <- lines[c(1,3,5,7)+1]
lines  <- lines[c(1,2)]
#line + group
#x <- x[x$line %in% lines, ]
#lines <- unique(x$)
#lines <- lines[c(1,3,5,7,9,11,13,15)]
#lines <- lines[c(1,3,5,7,9,11,13,15)+1]
list_A <- list()
list_B <- list()
for(i in 1:length(lines)){
  #i <- i
  #x_i <- colSums(x[which(x$line==lines[[i]]),-c(6222)])
  #x_i <- colSums(x[which(x$line==lines[[i]]),-c(953)])
  x_i <- colSums(x[which(x$line==lines[[i]]),-c(ncol(x))])
  x_c <- as.data.frame(x_i)
  colnames(x_c) <- c(lines[[i]])
  list_B[[i]] <- x_c
  x_i <- x_i[x_i >0]
  x_i <- names(x_i)
  list_A[[i]] <- x_i
}

list_B <- do.call(cbind,list_B)
require(RColorBrewer)
#heatmap(as.matrix(list_B),cexCol = 0.7,cexRow = 0.3,Rowv = NA,na.rm = T)
names(list_A) <- lines

lapply(1:length(list_A),function(i){
  print(paste(names(list_A)[[i]],length(list_A[[i]])))
}) 

ggVennDiagram::ggVennDiagram(list_A,set_size = 2,#3,#, 
                             #force_upset = F#, 
                             #order.set.by = "name", 
                             #order.intersect.by = "none"
                             ) +
  scale_fill_distiller(palette = "Reds", direction = 1)
#venn = ggVennDiagram::Venn(list_A)
#ggVennDiagram::plot_upset(venn, nintersects = 10, relative_height = 3, relative_width = 1,top.bar.numbers.size = 4,)
# lines <- factor(lines,lines[c(13,14,5,6, #l8
#                               15,16,7,8,#l9
#                               9,10,1,2,#l12
#                               11,12,3,4#l14
# )])
#list_A
#2021 not in 2022
table_Tax <- as.data.frame(tax_table(psd2))
#2021
X2021 <- table_Tax[(row.names(table_Tax) %in% list_A[[1]][!list_A[[1]] %in% list_A[[2]]]),]
X2021$year <- 2021
#2022
X2022 <- table_Tax[(row.names(table_Tax) %in% list_A[[2]][!list_A[[2]] %in% list_A[[1]]]),]
X2022$year <- 2022
##
x_dif <- rbind(X2021,X2022)
write.csv(x_dif,"C:/TESIS/DIF_HEAT.csv",row.names = F)
upset(fromList(list_A), order.by = "freq", mb.ratio = c(0.55, 0.45),nsets = 16,
      #sets=lines,
      decreasing = T,nintersects = 50,keep.order = TRUE)
################################################################################



otu_jac[metadata$SampleID[which(metadata$group==groups[[1]])],]

jacc_dist <- vegdist(otu_jac,method = "jaccard")

plot(hclust(jacc_dist))
# lines <- unique(metadata$group)
# list_AA <- list()
# for(i in 1:length(lines)){
#   
#   list_AA[[i]] <- colnames(otu_jac[metadata$SampleID[which(metadata$group==lines[[i]])],])
# 
# }
# names(list_AA) <- lines
# 
# venn = ggVennDiagram::Venn(list_AA)
# ggVennDiagram::plot_upset(venn, nintersects = 10, relative_height = 3, relative_width = 1,top.bar.numbers.size = 4,)



# 
# 
# 
# 
# 
# glom <- tax_glom(psd2, taxrank = 'Phylum', NArm = T)
# x <- as.data.frame(glom@otu_table)
# 
# x[x>0] <- 1
# colnames(x) <-glom@tax_table[,2]
# #colnames(x) <- glom@tax_table[,5]#5 family, 6 genus
# #x$line <- metadata$line
# #x$line <- paste0(metadata$line,"-",metadata$group)
# x$line <- glom@sam_data$group
# #x$line <- metadata$group
# #lines <- unique(x$line)
# 
# #line

otu_jac2 <- otu_jac
otu_jac2$line <- psd2@sam_data$group
tax_phyla <- cbind(row.names(psd2@tax_table),psd2@tax_table[,2])
phyla <- unique(tax_phyla[,2])
lines <- unique(psd2@sam_data$group)

list_AA <- list()
for(i in 1:length(lines)){
# i <- 1
  #i <- 1
  x_i <- otu_jac2[which(otu_jac2$line==lines[[i]]),]
  list_AA_j<- list()
    for(j in 1:length(phyla)){
      #j <- 1
      #j <- 23
      col_to <- tax_phyla[,1][which(tax_phyla[,2]==phyla[[j]])]
      #n sp
      if(length(col_to)==0){
        count <- 0
      } else if(length(col_to)==1){
        if(length(sum(x_i[,col_to]))==0){
          count <-0
        } else {
          count <- 1
        }
        count
      } else {
        count <- length(colSums(x_i[,col_to]) >0)  
      }
      
      list_AA_j[[j]] <- data.frame(phylum =phyla[[j]], 
                                 count=count)
    }
  list_AA_j <- do.call(rbind,list_AA_j)
  #list_AA_j$line <- lines[[i]]
  list_AA[[i]] <- list_AA_j 
  rm(list_AA_j)
}

list_AA <- do.call(cbind,list_AA)

list_AA <- list_AA[,-c(3,5,7)]
colnames(list_AA) <- c("Phylum",lines)
############################
############################
x_r <- as.data.frame(psd2r@otu_table)
x_r$line <- paste0(psd2r@sam_data$line,"-",psd2@sam_data$group)

lines <- unique(x_r$line)


list_lines <- list(
  line1 = c("L12-heat_2021","L12-control_2021"),
  line2 = c("L12-heat_2022","L12-control_2022"),
  line3 = c("L14-heat_2021","L14-control_2021"),
  line4 = c("L14-heat_2022","L14-control_2022"),
  line5 = c("L8-heat_2021","L8-control_2021"),
  line6 = c("L8-heat_2022","L8-control_2022"),
  line7 = c("L9-heat_2021","L9-control_2021"),
  line8 = c("L9-heat_2022","L9-control_2022")
)


list_lines_l2 <- list()
for(i in 1:length(list_lines)){
  #i <- 1
  x1 <- log2(colMeans(x_r[which(x_r$line==list_lines[[i]][1]),-ncol(x_r)])/
               colMeans(x_r[which(x_r$line==list_lines[[i]][2]),-ncol(x_r)]))
  x1 <- x1[which(!is.infinite(x1))]
  x1 <- x1[which(abs(x1)>2)]
  x1 <- as.data.frame(x1)
  x1$ASV <- row.names(x1)
  x1 <- x1[,c(2,1)]
  colnames(x1) <- c("ASV","count")
  x1$group <- paste0(list_lines[[i]],collapse = "/")
  list_lines_l2[[i]] <- x1
}

taxx <- as.data.frame(psd2r@tax_table)
taxx$ASV <- row.names(taxx)
list_lines_l2 <- do.call(rbind,list_lines_l2)
list_lines_l2 <- left_join(list_lines_l2,taxx)
write.csv(list_lines_l2,"C:/TESIS/LOG2FOLD.csv",na = "",row.names = T)
