##clustering with our MS1 features output...
## file output is matchedFeatures.txt from MaxQuant
library(tidyverse)
library(gplots)
#library(dplyr)
library(DESeq2)
library(dendextend)
library(data.table)
library(cluster)
library(propagate)
library(MineICA)
#library(fastICA)
library(kmed)
library(here)
library(impute)
library(pvclust)
library(ggdendro)
library(cutr) #devtools::install_github("moodymudskipper/cutr")
#library(ggridges)
###
## Remove peptides that are only sometimes identified...
## Data filtering function
## https://datascienceplus.com/proteomics-data-analysis-2-3-data-filtering-and-missing-value-imputation/
### functions
is.nan.data.frame <- function(x){
    do.call(cbind, lapply(x, is.nan))}


filter_valids = function(df, conditions, min_count, at_least_one = TRUE) {
    # df = data frame containing LOG2 data for filtering and organized by data type
    # conditions = a character vector dictating the grouping
    # min_count = a numeric vector of the same length as "conditions" indicating the minimum 
    #     number of valid values for each condition for retention
    # at_least_one = TRUE means to keep the row if min_count is met for at least one condition
    #     FALSE means min_count must be met across all conditions for retention
    #df[df==0] <- NA
    df[is.nan.data.frame(df)] <- 0
    all_names <- colnames(df) 
    cond.names = lapply(conditions, # Group column names by conditions
                        function(x) grep(x, all_names, value = TRUE, perl = TRUE))
    cond.filter = sapply(1:length(cond.names), function(i) {
        df2 = df %>% dplyr::select(cond.names[[i]])   # Extract columns of interest
        df2 = as.matrix(df2)   # Cast as matrix for the following command
        sums = rowSums(df2!=0) # count the number of valid values for each condition
        sums >= min_count[i]   # Calculates whether min_count requirement is met
    })
    if (at_least_one) {
        df$KEEP = apply(cond.filter, 1, any)
    } else {
        df$KEEP = apply(cond.filter, 1, all)
    }
    #return(df) # No rows are omitted, filter rules are listed in the KEEP column
    
    return(df)  # only keeping rows that meet the criteria!
}

filter_valids_quart = function(df, conditions, min_count, int_remove, at_least_one = TRUE) {
    # df = data frame containing LOG2 data for filtering and organized by data type
    # conditions = a character vector dictating the grouping
    # min_count = a numeric vector of the same length as "conditions" indicating the minimum 
    #     number of valid values for each condition for retention
    # at_least_one = TRUE means to keep the row if min_count is met for at least one condition
    #     FALSE means min_count must be met across all conditions for retention
    #df[df==0] <- NA
    df[is.nan.data.frame(df)] <- 0
    df[df <= int_remove] <- 0
    all_names <- colnames(df) 
    cond.names = lapply(conditions, # Group column names by conditions
                        function(x) grep(x, all_names, value = TRUE, perl = TRUE))
    cond.filter = sapply(1:length(cond.names), function(i) {
        df2 = df %>% dplyr::select(cond.names[[i]])   # Extract columns of interest
        df2 = as.matrix(df2)   # Cast as matrix for the following command
        sums = rowSums(df2!=0) # count the number of valid values for each condition
        sums >= min_count[i]   # Calculates whether min_count requirement is met
    })
    if (at_least_one) {
        df$KEEP = apply(cond.filter, 1, any)
    } else {
        df$KEEP = apply(cond.filter, 1, all)
    }
    #return(df) # No rows are omitted, filter rules are listed in the KEEP column
    
    return(df)  # only keeping rows that meet the criteria!
}




### For the openMS files
## will have to edit...
## sample info ../data/ms1_openms_sampleinfo.csv
ms1_openms <- fread(here::here("data", "ms1_openms_int.txt"), sep='\t')
ms1_openms_samples <- read.csv(here::here("data", "ms1_openms_sampleinfo.csv"), header=F)
ms1_openms_metadata <- read.csv(here::here("data", "plate2_drug_metadata_CPHSWAP.csv")) %>%
    mutate(plate_cell = paste0("P2", Cell))

#ms1_openms_metadata <- read.csv(here::here("data", "ms1_openms_metadata.csv")) %>%
#    mutate(plate_cell = paste0("P2", Cell))


#ms1_header <- read.table("./5.Zhibin_knime_output/combined_features_header.txt", header=F)
#ms1 <- read.table("./4.MS1features_linked/ms1features_consensus_minimal_edit.csv", sep=",",
#                 header=T) %>% as.data.frame() %>% dplyr::select(-CONSENSUS)

## Grabbing names of samples
#sample_names <- read.table("./5.Zhibin_knime_output/sample_names.txt", header=F) %>% 
#    mutate(Sample_name = substr(V1, 1, nchar(V1 %>% as.character())-5))
#map_names <- read.table("./5.Zhibin_knime_output/map_names.csv", sep=",", header=T) 
#map_names <- map_names %>% mutate(Sample_name = substr(File, 1, nchar(File %>% as.character())-13)) %>%
#    merge(., sample_names, by="Sample_name", all.y=T)


consensus_features <- ms1_openms %>% dplyr::select(.,contains("cf"))
ms1 <- ms1_openms %>% dplyr::select(.,-contains("cf"))
rt <- ms1 %>% dplyr::select(.,contains("rt"))
mz <- ms1 %>% dplyr::select(.,contains("mz"))
charge <-  ms1 %>% dplyr::select(.,contains("charge"))
intensity <- ms1 %>% dplyr::select(., contains("intensity")) %>% as.data.frame()
colnames(intensity) <- ms1_openms_samples$V1

#write.csv(intensity, file="~/iclouddrive/Documents/Manuscripts/figeys/clusterMethods/manuscript/clust_rename/data/ms1/ms1_intensity.csv")


int_col <- colnames(intensity) %>% data.frame(.)

ms1_openms_metadata <- left_join(int_col, ms1_openms_metadata, by=c("." = "plate_cell")) ## only including samples we used
## try to reduce number of features...
# reorganize df so that it's in the order of drugs (for filtering)
meta_drug <- ms1_openms_metadata %>% arrange(., Drug)
intensity <- intensity %>% dplyr::select(meta_drug$.)
colnames(intensity) <- meta_drug$Name
#intensitydf <- as.data.frame(intensity)
#intensitydf <- reshape2::melt(intensitydf, value.name = "intensity") # reshape the matrix 

#ggplot(intensitydf %>% filter(intensity != 0), aes(x = intensity)) +
#    geom_density(alpha = 0.2, size = 1.25) 

cond_options_x <- table(ms1_openms_metadata$Drug) %>% as.data.frame()

## options for drug specific filtering
## want to change intensities below a noise threshold for 
cond_opts <- cond_options_x$Var1
cond_count <- cond_options_x$Freq * 0.5

## lets look at the distribution of intensities, let's consider the bottom quartile...so if intensity is in bottom quartile, 
## we will consider if as a "0" ONLY in our filtering step!

int <- intensity %>% pivot_longer(names_to = "Name", values_to="intensity", cols=everything()) %>% drop_na()
summary(int$intensity)
int_4_quartiles <- smart_cut(int$intensity,4,"groups", labels = c("remove", "Low", "Med", "High"))
table(int_4_quartiles)
int$quart <- int_4_quartiles
max_to_remove <- int %>% filter(quart == "remove")
max_to_remove <- max_to_remove$intensity %>% max()


### Can we visualize the distribution of each intensity?
head(int) ## want to break it up by quartiles

#ggplot(data = int,
#       aes(x = log2(intensity), fill = quart)) +
#    geom_density_ridges() + #(alpha = .7, rel_min_height = .01,
#                        #color = "white", from = -5, to = 95) +
#    theme_ridges(grid = FALSE) 
#    
#ggplot() +
#    geom_density(data=int, aes(x=log2(intensity), fill=as.factor(quart)),
#                 colour="black", alpha=0.6)  +
#    xlab("Peptide Intensity") +
#    ylab("Density") +
#    theme(legend.position="bottom", legend.key = element_rect(colour = "black")) 
#    #scale_x_continuous(breaks=seq(0,300, by=30), limits=c(0,300)) 



dens <- density(log2(int$intensity))
df <- data.frame(x=dens$x, y=dens$y)
probs <- c(0,0.25,0.5,0.75,1)
quantiles <- quantile(log2(int$intensity), prob=probs)
df$quant <- factor(findInterval(df$x,quantiles))


(pep_quart_plot <- ggplot(df, aes(x,y)) + geom_line() + 
    geom_ribbon(aes(ymin=0, ymax=y, fill=quant)) + 
    scale_x_continuous(breaks=quantiles) + 
    #scale_fill_brewer(guide="none") +
    xlab("Peptide intensity quartiles, log2") +
    ylab("Density") +
    scale_fill_viridis_d(option = "plasma", guide = "none") +
    theme(text=element_text(size=12))
)
#ggsave(here::here("output", "figs", "ms1_pretest_quart_swap_intquart.pdf"), pep_quart_plot, 
 #      width = 8, height = 4, units ="in")




## filter by at least Q50 in the treatments
drug_exp_filt <- filter_valids(intensity, 
                               conditions = cond_opts,
                               min_count = cond_count, 
                               at_least_one = T)  ## THERE ARE NO MISSING VALUES!!!!
drug_exp_filt <- intensity[drug_exp_filt$KEEP == 1,] ## 48,543

## filter by at least Q50 in each treatment considering ZEROS and anything in low quartile
drug_exp_filt_quart <- filter_valids_quart(intensity, 
                               conditions = cond_opts,
                               min_count = cond_count, 
                               int_remove = max_to_remove,
                               at_least_one = T)  ## THERE ARE NO MISSING VALUES!!!!
drug_exp_filt_quart <- intensity[drug_exp_filt_quart$KEEP == 1,] ## 37,484

## Normalize by library size
norm_pep <- estimateSizeFactorsForMatrix(drug_exp_filt_quart)
norm_exp <- sweep(drug_exp_filt_quart, 2, norm_pep, "/") ##peptides are normalized

## impute missing values
norm_exp[norm_exp== 0] <- NaN
#norm_exp_noCPH[norm_exp_noCPH==0] <- NA

imputed_pep<- impute.knn(norm_exp %>% as.matrix(), k = 10, rowmax = 0.5, colmax = 0.95, rng.seed = 362436069)


#norm_exp_noCPH <- norm_exp %>% dplyr::select(-contains("CP-H"))

### uncomment to write the filtered data to file
#write.csv(norm_exp, file=here::here("data","ms1_openms_filt_quart.csv"))

## do the 
#write.csv(norm_exp_noCPH, file=here::here("data","ms1_openms_filt_quart-NO-CP.csv"))


## Read in KNN from file
## from python knn
#ms1_knn  <- fread(here::here("data", "ms1_openms_quart.knn.csv"), sep=",", header = F)
#ms1_knn  <- fread(here::here("data", "ms1_openms_quart-NO-CP-knn.csv"), sep=",", header = F)
ms1_knn <- imputed_pep$data

#colnames(ms1_knn) <- colnames(drug_exp_filt_quart)
colnames(ms1_knn) <- colnames(norm_exp)




## impute missing values
#norm_exp[norm_exp== 0] <- NaN
#imputed_ms1<- impute.knn(norm_exp %>% as.matrix(), k = 10, rowmax = 0.5, colmax = 0.9, rng.seed = 362436069)
#


## NO NORMALIZATION!
## should we normalize differently?
#norm_int <- estimateSizeFactorsForMatrix(ms1_knn) 
#norm_data <- sweep(ms1_knn, 2, norm_int, "/") ##peptides are normalized

# Transform by log2
ms1_log <- data.frame(ms1_knn) %>% #dplyr::select(starts_with('Intensity')) %>%
    mutate_all(., funs(log2(.)))
colnames(ms1_log) = colnames(ms1_knn)
#log_exp[log_exp == 0 ] <- 0.001 # don't do this for changing zeros



ctrl <- dplyr::select(ms1_log, contains("NC"))
ctrl$median = apply(ctrl, 1, median, na.rm = T)
drug_temp <- dplyr::select(ms1_log, -contains("NC"))
ctrl<- matrix(ctrl$median, nrow=length(ctrl$median), ncol=ncol(drug_temp), byrow=F)


fc <- ctrl - drug_temp
fc$feat_num <- rownames(drug_exp_filt_quart)

## want to use as many components in ICA as what explains 99% of variation in PCA...
noctrlmeta <- ms1_openms_metadata %>% filter(Drug != "NC") %>%
    arrange(., Drug)
#noctrlCPHmeta <- ms1_openms_metadata %>% filter(Drug != "NC" & Meta != "CP-H-48")

fc <- fc %>% dplyr::select(noctrlmeta$Name)
#fc <- fc %>% dplyr::select(noctrlCPHmeta$Name)

#write.csv(fc, here::here("output", "ms1_openms_quart", "ms1_openms_quart_log2fc.txt"), quote=F, row.names = T)

#write.csv(fc, here::here("output", "ms1_openms_quart_noCPH", "ms1_openms_quart_noCPH_log2fc.txt"), quote=F, row.names = T)
#write.csv(fc, here::here("output", "ms1_openms_quart_swap", "ms1_openms_quart_swap_log2fc.txt"), quote=F, row.names = T)
pca<- prcomp(fc %>% t(), center=T, scale=F)
sampleVals<-data.frame(pca$x)
exprVals<-data.frame(pca$rotation)
PoV <- (pca$sdev^2/sum(pca$sdev^2))*100 



coords<-data.frame(sampleVals, Drug = noctrlmeta$Drug,
                   samplename = rownames(sampleVals),
                   Concentration = noctrlmeta$Conc)
coords$Concentration <- factor(coords$Concentration, levels = c("H", "M", "L"))
numPCs <- 1:length(PoV)

for (i in 1:length(PoV)) {
    percent <- paste0("(", round(PoV[i],2), "%)")
    name <- paste0("PC", i, "per")
    assign(name, percent)
}

quick_pca <- function(coords_df, x_PC, y_PC, shape, fill){
    yaxis <- y_PC
    xaxis <- x_PC
    yperc <- paste0("(", round(PoV[yaxis %>% as.numeric()] ,2), "%)")
    xperc <- paste0("(", round(PoV[xaxis %>% as.numeric()] ,2), "%)")
    
    yaxislabel <- paste0("PC", yaxis, " ", yperc)
    xaxislabel <- paste0("PC", xaxis, " ", xperc)
    pca_plot <- ggplot(coords, aes_string(x = paste0('PC', xaxis), y = paste0('PC', yaxis))) + #accept selectInput to choose axes!
        geom_point(size=4, aes_string(shape = shape, color = fill)) + 
        #scale_color_manual(values=c("#3444DA",  "#FFBA49", "#8C5999", "#D7C0D0", "#363635", "#DF737A", "#20A39E", "#8EB8E5")) + 
        scale_color_manual(values=c("#3444DA",  "#DF737A", "#8C5999", "#D7C0D0", "#FFBA49")) +
        scale_x_continuous(name=xaxislabel) + # labels depend on selected PCs
        scale_y_continuous(name=yaxislabel) +
        theme_bw() +
        theme(legend.position = "right", text=element_text(size=16))
    return(pca_plot)
}
(pc12 <- quick_pca(coords_df = coords, x_PC = 1, y_PC = 2, shape = "Concentration", fill = "Drug"))

(pc23 <- quick_pca(coords_df = coords, x_PC = 2, y_PC = 3, shape = "Concentration", fill = "Drug"))


#ggsave(here::here("output", "figs", "ms1_pretest_quart_swap_pc12.pdf"), pc12, width = 7, height = 5, units = "in")
#ggsave(here::here("output", "figs", "ms1_pretest_quart_swap_pc23.pdf"), pc23, width = 7, height = 5, units = "in")

## Want to cluster the logFC for each individual by DMSO



## consider the clusters 

################################
## ICA and kmedoids in python ##
################################
sil <- read.csv(here::here("output", "ms1_openms_quart_swap", "k-med_silhouette.csv"), header=F)
inertia <- read.csv(here::here("output","ms1_openms_quart_swap", "k-med_inertia.csv"), header=F)
kmedoid_labels <- read.csv(here::here("output", "ms1_openms_quart_swap", "k-med_labels.csv"), header=F) %>% as.vector()

#rownames(kmedoid_labels) <- filt_peptides
clusterinfopep <- kmedoid_labels %>% data.frame()

max_clust <- kmedoid_labels %>% unique() %>% max()
fc_clusters <- list()
for (x in 0:max_clust){
    clustername <- paste("Cluster", x)
    y <- kmedoid_labels[which(kmedoid_labels$V1 == x),,drop=F] %>% rownames()
    fc_clusters[[clustername]] <- fc[y %>% as.numeric(),]
    fc_clusters[[clustername]]$feat_num <- y 
}

### eigengenes

### Singular value decomposition for genome-wide expression data processing and modeling
###  https://doi.org/10.1073/pnas.97.18.10101
## https://www.biostars.org/p/299167/
## eigengenes <- svd(X)$v
### but WGCNA calculates the first principal component...so whats the difference?
pca_clusters <- lapply(fc_clusters, function(x) prcomp(x %>% dplyr::select(-feat_num) %>% t(), center=T, scale=F)$x)
pc1_cluster <- lapply(pca_clusters,  "[", , 1)

clustereigengenes <- bind_cols(pc1_cluster) %>% data.frame()
### make a factored dataframe of the compounds
## make a factored dataframe of the compounds

noctrlmeta <- noctrlmeta %>% mutate(drug_conc = paste0(Drug, "-", Conc))

datTraits <- data.frame(compounds= as.factor(noctrlmeta$drug_conc))
rownames(datTraits) <- noctrlmeta$Name
traits <- model.matrix(~ ., data=datTraits,
                       contrasts.arg = lapply(datTraits, contrasts, contrasts=FALSE))
traits <- traits[,-1] ##removing intercept
moduleTraitCor <- stats::cor(clustereigengenes, traits, method = "s")
#moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,nrow(norm_data))
#
colnames(moduleTraitCor) <- colnames(moduleTraitCor) %>% substr(., 10, nchar(colnames(moduleTraitCor)))

###hclust clusters the rows
module_clust <- hclust(as.dist(1-cor(t(moduleTraitCor), method="pearson")), method="average")
drug_clust <- hclust(as.dist(1-cor(moduleTraitCor, method="pearson")), method="average")
drug_clust_boot <- pvclust(moduleTraitCor, method.hclust="average",
                           method.dist="correlation", nboot=1000,
                           iseed=2248)

#plot(drug_clust)
#pvrect(drug_clust_boot, alpha=0.9, pv="au")
clust_boot <- as.dendrogram(drug_clust_boot)
drug_clust_boot %>% as.dendrogram() %>% plot(main = "test")
drug_clust_boot %>% text 
drug_clust_boot %>% pvrect(alpha=0.9, pv="au")


# colour palette for heatmap
my_palette <- colorRampPalette(c("#8C5999", "white", "#FEAD34"))(n = 299)  #FF7620 #blue is 3444DA 
mango_balm <- colorRampPalette(c("#6D8CB9", "white", "#FDB62A"))(n=200)
new_palette <- colorRampPalette(c("#70A9A1", "white", "#FF9633"))(n=120)
dark_palette <- colorRampPalette(c("#70A9A1", "white", "#FF9633"))(n=299)
# colour for treatment info beneath plot
## "#3444DA",  "#DF737A", "#8C5999", "#D7C0D0", "#FFBA49"
cor_meta <- data.frame(treatment = colnames(moduleTraitCor)) %>%
    mutate(drug_col = case_when(str_detect(treatment, 'AZ') ~ '#3444DA',
                                str_detect(treatment, 'CP') ~ '#DF737A',
                                str_detect(treatment, 'DC') ~ '#8C5999',
                                str_detect(treatment, 'NZ') ~ '#D7C0D0',
                                str_detect(treatment, 'PR') ~ '#FFBA49'
    )) %>%
    mutate(conc_col = case_when(str_detect(treatment, '-H') ~ '#08090A',
                                str_detect(treatment, '-M') ~ '#96999C',
                                str_detect(treatment, '-L') ~ '#EAEBEB'))
#cor_meta <- cor_meta[drug_clust$order,]

#pdf(here::here("output","figs", "pretest-ms1_fc_kmed_hm2-quar_swap.pdf"), width=9, height=6)

# if getting "invalid graphics state" error, try:
# par(mar=c(1,1,1,1))
par(xpd = TRUE) # allows lgend to be outside "official" plotting lines
coords2 <- list(x=0, y=0.95)
heatmap.2(moduleTraitCor,
          #cellnote = corr_clust,  # same data set for cell labels
          #main = "Correlation", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,0),     # widens margins around plot
          #col=my_palette,       # use on color palette defined earlier
          #col=mango_balm,
          col=new_palette,
          breaks=seq(-0.6,0.6,0.01),    # enable color transition at specified limits
          dendrogram="col",     # only draw a row dendrogram
          Colv=as.dendrogram(drug_clust),            # turn off column clustering
          #Colv = test_dend,
          Rowv=as.dendrogram(module_clust),
          #ColSideColors = cor_meta$drug_col,
          ColSideColors = cor_meta$drug_col,
          labRow = "",
          key.xlab ="Correlation",
          key.title = NA,
          lhei = c(1,4),
          lwid= c(1,3),
          cexCol = 1.2,
          keysize = 1)
#coords=locator(1)
legend(coords2, title = "Drug",legend=c("AZ","CP", "DC", "NZ", "PR"), 
       fill=c("#3444DA",  "#DF737A", "#8C5999", "#D7C0D0", "#FFBA49"), cex=0.8, box.lty=0)
#dev.off()


## save RData image for dendrogram comparison
ms1_moduleTraitCor <- moduleTraitCor
ms1_drugclustboot <- drug_clust_boot
ms1_drugclust <- drug_clust 
#save(ms1_moduleTraitCor, ms1_drugclust, ms1_drugclustboot, file=here::here("output", "ms1_openms_quart_swap", "ms1_dedro.RData"))

#
###### This will be important later on! Only consider significant correlations
##moduleTraitCor[moduleTraitPvalue > 0.05] <- NA
#
#### reorder based on clustering
#corr_clust <- moduleTraitCor[module_clust$order, drug_clust$order]
#
##
###GGplot style heatmap
#correlation_ggplot <- corr_clust%>% reshape2::melt()
##
#dx <- dendro_data(drug_clust)
#dy <- dendro_data(module_clust)
###
#
### helper function for creating dendograms
#ggdend <- function(df) {
#    ggplot() +
#        geom_segment(data = df, aes(x=x, y=y, xend=xend, yend=yend)) +
#        labs(x = "", y = "") + theme_minimal() +
#        theme(axis.text = element_blank(), axis.ticks = element_blank(),
#              panel.grid = element_blank())
#}
##
### x/y dendograms
###library(dendextend)
#px <- ggdend(dx$segments)
#py <- ggdend(dy$segments) + coord_flip()
#
#
### convert this into pheatmap
## https://www.datanovia.com/en/lessons/heatmap-in-r-static-and-interactive-visualization/
#(heatmap_full <- ggplot(data = correlation_ggplot, 
#                        aes(x = Var2, y = Var1)) +
#        geom_tile(aes(fill=value), colour="white") +
#        scale_fill_gradientn(colours=c("#344499", "#59C3C3","white","#DF737B",  "#533747"),
#                             space = "Lab", name="Correlation", na.value = "grey95",
#                             limits=c(-.45,0.45)) +
#        #coord_flip() +
#        labs(x="Drug treatment", y="Clusters") +  theme_bw() + #theme(axis.text.y=element_blank(), plot.margin = margin(-0.75, 0, 0,0 , "cm")))
#        theme(axis.text.y=element_text(angle=20,vjust=0.5),
#              axis.text.x=element_text(angle=90, vjust=0.5),
#              plot.margin = margin(0, 0, 0,0 , "cm")))
#ggsave(here::here("output","figs", "ms1openms_fc_kmed_hm.pdf"),heatmap_full, width=8, height=7, units="in" )
#ggsave(here::here("output","figs", "ms1openms_fc_kmed_dendro.pdf"), px, width=8, height=1,  units="in" )
#