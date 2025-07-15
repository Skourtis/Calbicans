# check bicor against GS
library(data.table);library(ggplot2);library(WGCNA)
library(PRROC);library(ggpubr);library(stringr)
# Uniprot mapping
Uniprot_annot = fread(here::here('in','datasets','uniprotkb_taxonomy_id_237561_2025_07_08.tsv.gz'))
Uniprot_annot[,CGD:=str_remove(CGD,';[:print:]*$')]

# check the covariation strength of different Go terms/
# subcellular compartments across subsets
# similar to complexes

# GoldStandard
GS = fread(here::here('out','datasets','STRING_TP_FP.gz') )
#### Annotate with gold standard ####
setnames(GS,c('Protein_1','Protein_2','Class'))

  # isolates with different stress
  subset_stress = fread(here::here('in','datasets','subset_conditions_proteomes_refreplicates.tsv'))
  subset_stress = subset_stress[!is.na(Genes) & plate_position != 'qc']
  subset_stress[,plate_position:=paste(plate_position,condition,sep = '_')]
  # why some have the same names
  subset_stress = subset_stress[,.(PG.MaxLFQ = mean(PG.MaxLFQ,na.rm = T)),by= .(plate_position,Genes)]
  subset_stress[,.(.N),by = .(Genes,plate_position)][,N] |> hist()
  subset_stress[,abundance:= log2(PG.MaxLFQ)]
  subset_stress[,mean_expression:=median(abundance,na.rm = T),by = .(plate_position)]
  subset_stress[,abundance := abundance -mean_expression ]
  
  subset_stress = subset_stress|> dcast(Genes~plate_position,value.var = 'abundance') 
  subset_stress |> nrow()
  subset_stress = merge(subset_stress,Uniprot_annot[,.(Entry,CGD)],by.x = 'Genes',by.y = 'CGD')
  subset_stress |> nrow() 
  subset_stress[,Genes:=NULL]
  subset_stress <- subset_stress |> tibble::column_to_rownames('Entry')
  subset_stress |> dim()
  df = subset_stress

  
  
  # Transpose and median normalise the rows to avoid spurious correlations
  tmp_medians <- apply( df, 1, median, na.rm = TRUE )  
  df <- sweep( df, 1, tmp_medians, FUN = "-" )
  
  
  pca_res <- prcomp(df  %>% na.omit() %>% t(), scale=F)
  var_explained <- pca_res$sdev^2/sum(pca_res$sdev^2)
  
  df_tmp = pca_res$x |> as.data.table()
  df_tmp$sample = rownames(pca_res$x) 
  df_tmp$condition = df_tmp$sample  |> str_extract('_S[:print:]*') |> str_remove('^_')
  df_tmp[,N_samples := .N, by = condition]
  pca_plot =   ggplot(df_tmp,aes(x=PC1,y=PC2, colour = condition)) + 
    geom_point(size=1)+
    theme_bw()+
    # geom_point()+
    # geom_point(data = df_tmp[condition %in% c('PXD024229','PXD007705','PXD006383','PXD013496','PXD002065')], 
    #            aes(colour = project))+
    theme(legend.position = 'bottom')
  
  ggsave(here::here('out','plots','pca_samples_plot.pdf'),pca_plot)
  
  list_dataset = list()
  for(i in unique(df_tmp$condition )){
    print(i)
    
    list_dataset[[i]]= df[,str_detect(colnames(df),paste0(i,'$'))]
    
  }
  list_dataset |> purrr::map_dbl(ncol)
  
  
list_PRcurves = list()
n_detected = 30
for(dataset in names(list_dataset)){

  df_subset = list_dataset[[dataset]] 
#### Calculate correlations between proteins ####
# Select only proteins with at least ratios
feature_count <- apply(df_subset, 1, function(x){ sum(!is.na(x))} )
# almost all proteins detected in all samples
feature_count |> hist()
df_subset <- df_subset[feature_count >= n_detected,]
# Transpose and median normalise the rows to avoid spurious correlations
tmp_medians <- apply( df_subset, 1, median, na.rm = TRUE )  
df_subset <- sweep( df_subset, 1, tmp_medians, FUN = "-" )
df_subset = df_subset|> t()
# Define function for bicor correlation. It takes an input data frame and a name of the output value. 
# It outputs a row-wise sorted list of protein pairs along with calculated distance metric
f_BIC <- function( input_data, value_name ){
  missingness_tmp <- base::crossprod(!is.na(input_data))>30 # with microproteins being very rare it's common to only have 1 point in common produce fake correlations
  #using suggested bicor
  tmp <- bicor( input_data , use = "pairwise.complete.obs", nThreads = 6)                     # Robust correlation using all pairwise complete observations
  tmp[missingness_tmp == FALSE] <- NaN
  tmp <- as.data.table( reshape2::melt( as.matrix( tmp  )))                                   # Turn distance matrix into a pair-wise data.table
  tmp <- tmp[, .( Protein_1 = as.character(Var1), Protein_2 = as.character(Var2), value ) ]   # Rename and change to character
  tmp <- tmp[ Protein_1 > Protein_2 ]                                                         # Remove redundant pairs by keeping only A > B, removing A == B and B < A pairs
  names(tmp)[3] <- value_name                                                                 # Assign new value name
  return(tmp)
}

# Calculate a set of pairwise similarities
list_PRcurves[[dataset]] <- f_BIC(df_subset, dataset)

}
# Combine the results
DT <- Reduce(function(...) merge(...,by = c("Protein_1", "Protein_2")),list_PRcurves)    # This works because all pairs are sorted A > B
DT$avg_covariation = DT[,-c(1:2)] |> as.matrix() |> matrixStats::rowMeans2(na.rm = T)
# DT$KO |> hist()

# Simplify protein IDs for merging with gold standard
DT <- DT[ Protein_1 != Protein_2 ]                  # Removing self-associations that arose between isoforms (can't test those anyway)
DT[ Protein_1 <= Protein_2 ]                        # Validate that A > B still fullfilled   

# Annotate correlations with gold standards by merging
DT <- merge(DT, GS, by = c("Protein_1", "Protein_2"))
DT[, .N, Class]  # Currently about 3% TP pairs

# Downsample FPs to 90% - this is essential to make curves comparable across project
n_FP_to_remove <- DT[ Class == "FP", .N ] - DT[ Class == "TP", .N ] * 9 
DT[ sample( which( Class == "FP" ), n_FP_to_remove ), Class := NA ]    # From the rows where the class is FP, randomly sample a certain subset and turn to NA
DT <- DT[ !is.na(Class)  ][Class != '']                                             # Remove those rows
DT[, .N, Class]  # Now 10% TP pairs

# Add a randomised classifier to the input table
DT[, random_classifier := sample(`SC`)]


#### Run precision recall analyses ####

#### Define PR curve function and PR input data sets ####

# An R function for precision - recall analysis using the PRROC package.
# Inputs: data = a data table to analyse; gs = the column containing the class labels; metric = the column containing the metric to be tested
# Optional input: the "modifier" argument specifies how pairs will be ranked. The default "descending" means that pairs are sorted such that higher values correspond to higher confidence of association. For other metrics (e.g. distances), set to "ascending" to reverse the order, or to "absolute" to rank pairs by the absolute metric value
# Output: A list of length 2, containing some stats and the actual curve coordinates, respectively
PR_curve <- function(data, gs, metric, modifier = "descending"){
  tmp <- data[ !is.na(get(gs)) & !is.na(get(metric)), .(tmp_class = get(gs), tmp_metric = get(metric)) ]   # Subset the data to remove NAs, select and re-name the required columns
  if( modifier == "ascending" ){ tmp[, tmp_metric := -tmp_metric] }          # Reverse the tmp_metric if modifier set to "ascending"
  if( modifier == "absolute" ){ tmp[, tmp_metric := abs(tmp_metric)] }       # Obtains absolute values of tmp_metric if modifier set to "absolute"
  tmp_TP <- tmp[ tmp_class == "TP" , tmp_metric ]                            # Scores for TP pairs
  tmp_FP <- tmp[ tmp_class == "FP" , tmp_metric ]                            # Scores for TP pairs
  # tmp_FP <- sample( tmp_FP , length(tmp_TP)*9 )                              # Down-sample FP pairs to 90%
  tmp_curve <- pr.curve( tmp_TP, tmp_FP, curve = TRUE, rand.compute = TRUE ) # Calculate PR curve details
  tmp_curve_dt <- data.table( tmp_curve$curve, gold_standard = gs, metric = metric )               # Extract curve coordinates and turn them into data.table
  setnames(tmp_curve_dt, old = c("V1", "V2", "V3"), new = c("Recall", "Precision", "Thresholds"))  # Add informative col names
  tmp_curve_dt[, Thresholds := NULL ]                                        # Don't need the thresholds
  tmp_curve_dt <- tmp_curve_dt[ Recall > 0.005 ]                             # Drop low recall points because they are pretty much random
  tmp_curve_dt <- rbind( tmp_curve_dt[ Recall  < 0.05                ][ sample(.N, 180) ],   # Downsample the curves for faster plotting. Keep more points for the low recall areas because they are less robust.
                         tmp_curve_dt[ Recall >= 0.05 & Recall < 0.1 ][ sample(.N, 200) ],
                         tmp_curve_dt[ Recall >= 0.10 & Recall < 0.2 ][ sample(.N, 200) ],
                         tmp_curve_dt[ Recall >= 0.20 & Recall < 0.4 ][ sample(.N, 200) ],
                         tmp_curve_dt[ Recall >= 0.40                ][ sample(.N, 200) ])
  output <- list( 
    "stats" = data.table( "AUPRC" = tmp_curve$auc.integral, "AUPRC_random" = tmp_curve$rand$auc.integral, "N_TP" = length(tmp_TP), gold_standard = gs, metric = metric), 
    "curve" = tmp_curve_dt)
  return(output)
}

# Perform the PR analyses using my custom function. 
PR_curves = list()
for(name in c(names(list_PRcurves),'avg_covariation','random_classifier')){
  print(name)
  PR_curves[[name]] = PR_curve( DT, "Class", name)
  
}

# Extract and combine PR statistics, then print them on screen
PR_stats = list()
PR_data = list()
for(name in c(names(list_PRcurves),'avg_covariation','random_classifier')){
  PR_stats[[name]] <- PR_curves[[name]][['stats']]
  PR_data[[name]] <- PR_curves[[name]][['curve']]
}

PR_stats <- Reduce(rbind,PR_stats)    # This works because all pairs are sorted A > B
PR_data <- Reduce(rbind,PR_data)    # This works because all pairs are sorted A > B

#### Plot the PR curves ####

# # Set standard plotting parameters
# plot_cols <- c( 
#   # prothd2_prohd1_subset = "#37abc8ff", 
#                 KO = "#ff2a7fff", random_classifier = "grey50") # Set plotting colours
# plot_labels <- c(
#   # prothd2_prohd1_subset = "proHD2_subset", 
#                  KO = "KO", random_classifier = "Random")
# 
# Define overall plotting theme
# theme_set( theme( line = element_line( linewidth = 0.25 ),
#                   panel.border = element_rect(colour = "black", linewidth = 0.25, fill = NA), panel.background = element_blank(), panel.grid = element_blank(),
#                   axis.text = element_text(size = 6), axis.title = element_text(size = 7),
#                   strip.background = element_blank(), strip.text = element_text(size = 8), 
#                   legend.title = element_blank(), 
#                   legend.background = element_blank(), legend.box.background = element_blank(), legend.text = element_text(size = 8), legend.key.width = unit(6, "mm"), 
#                   legend.key = element_rect(fill = NA, colour = NA),
#                   plot.title = element_text( size = 7, hjust = 0.5, vjust = 2)))
# 
# Define a function to create a PR curve plot with the AUPRC s8own as barchart as an inset
PR_plot <- function(PR_stats, PR_data, plot_title){
  
  # Plot the AUPRC bar-chart inset as a grob for custom_annotation
  inset <- ggplotGrob(ggplot(PR_stats, aes(x = metric, y = AUPRC, fill = metric))+
                        geom_bar(stat = "identity")+
                        scale_y_continuous( breaks = seq(0,1,0.2))+
                        # scale_fill_manual( values = plot_cols )+
                        theme( legend.position = "none", axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(),
                               axis.title.y = element_text(size = 8, angle = 90, vjust = 1), axis.text.y = element_text(size = 8)))
  
  # Create the actual plot
  p <- ggplot(PR_data, aes(x = Recall, y = Precision, colour = metric))+
    geom_line( linewidth = 0.5 )+
    ggtitle(plot_title)+
    # scale_colour_manual( values = plot_cols, labels = plot_labels )+
    scale_x_continuous( limits = c(0,1), breaks = seq(0,1,0.2))+
    scale_y_continuous( limits = c(0,1), breaks = seq(0,1,0.2))+
    # theme(legend.position = "none")+
    annotation_custom(grob = inset, xmin = 0.55, xmax = 1, ymin = 0.55, ymax = 1)
  
  return(p)
}
p <- PR_plot(PR_stats = PR_stats, PR_data = PR_data, 
             plot_title = glue::glue('Min_detection {n_detected} proteins'))+
  theme_bw()
ggsave(here::here('out','plots','STRING_GS_bicor_stressed_Subsets_PRcurves.pdf'),p)

# fit protein pairs which are changing
DT <- Reduce(function(...) merge(...,by = c("Protein_1", "Protein_2")),list_PRcurves)    # This works because all pairs are sorted A > B

DT = DT |> melt(id.vars = c('Protein_1','Protein_2'),variable.name = 'Condition',value.name = 'bicor')
DT[,.SD[sample(.N,10000)],by = .(Condition)] |> 
  ggplot(aes(y = bicor, x = Condition))+
  geom_boxplot()

# fit protein pairs which are changing
DT_normalised <- Reduce(function(...) merge(...,by = c("Protein_1", "Protein_2")),list_PRcurves)    # This works because all pairs are sorted A > B
colnames_samples = colnames(DT_normalised)
DT_normalised = cbind(DT_normalised[,1:2],
                      broman::normalize(DT_normalised[,-c(1:2)]))
colnames(DT_normalised) <- colnames_samples
DT_normalised = DT_normalised |> melt(id.vars = c('Protein_1','Protein_2'),variable.name = 'Condition',value.name = 'bicor')
DT_normalised[,.SD[sample(.N,10000)],by = .(Condition)] |> 
  ggplot(aes(y = bicor, x = Condition))+
  geom_boxplot()
zscore <- function(x) {
  (x - mean(x, na.rm=T) ) / sd(x, na.rm=T)
} 
DT_normalised[,`:=`(avg_cor = mean(bicor,na.rm = T),
                    sd_cor = sd(bicor,na.rm = T)),by = .(Protein_1,Protein_2)]
DT_normalised[,`:=`(condition_specific_cor = bicor-avg_cor)]
DT_normalised[,bands:=cut(avg_cor,1000)]
DT_normalised[,N_pairs_bands := .N, by = bands]
DT_normalised[, zscore := zscore(condition_specific_cor), by=bands]

fwrite(DT_normalised,here::here('out','datasets','covariation_stressed_treatments_normalised.gz'))
DT_normalised = fread(here::here('out','datasets','covariation_stressed_treatments_normalised.gz'))
across_cond_cor = DT_normalised |> 
  dcast(Protein_1 + Protein_2 ~ Condition, value.var = 'bicor') 
across_cond_cor |> dim()
across_cond_cor = f_BIC(across_cond_cor[,-c(1:2)],value_name = 'bicor')

across_cond_cor_rev = copy(across_cond_cor)
across_cond_cor_rev[,`:=`(Protein_1=Protein_2,
                      Protein_2=Protein_1)]
across_cond_cor = rbind(across_cond_cor,
                    across_cond_cor_rev)
across_cond_cor_matrix = dcast(across_cond_cor,Protein_1 ~ Protein_2, value.var = 'bicor') |> 
  tibble::column_to_rownames('Protein_1') |> 
  as.matrix()
diag(across_cond_cor_matrix) <- 1
pheatmap::pheatmap(across_cond_cor_matrix, main ='all protein pairwise covariations' )

# only for significant 
significant_pairs = merge(DT_normalised,
      unique(DT_normalised[avg_cor>0.57,.(Protein_1,Protein_2)])[,type:='significant'],
      all.x = T, by = c('Protein_1','Protein_2'))
significant_pairs[is.na(type), type:= 'non_significant']
significant_pairs[,head(.SD,1000000),by = .(Condition,type)] |> 
  ggplot(aes(x = bicor, fill =  type))+
  geom_density(alpha = 0.5)+
  facet_wrap('Condition')

library(caret)
across_cond_cor = DT_normalised |> 
  dcast(Protein_1 + Protein_2 ~ Condition, value.var = 'bicor') 
across_cond_cor |> dim()
across_cond_cor = across_cond_cor[,-c(1:2)]
across_cond_cor[across_cond_cor>0.95] <- T
across_cond_cor[across_cond_cor<=0.95] <- F
across_cond_cor = sapply(across_cond_cor, as.logical)
accuracy_matrix = matrix(rep(NaN,ncol(across_cond_cor)^2),
                         ncol = ncol(across_cond_cor))
for(i in 1:ncol(across_cond_cor)){
  print(i)
  for(j in 1:ncol(across_cond_cor)){
    accuracy = confusionMatrix(as.factor(across_cond_cor[,i]),
                               as.factor(across_cond_cor[,j]), 
                               positive = 'TRUE')
    accuracy_matrix[i,j] <- accuracy$overall[1]
  }
  
}
colnames(accuracy_matrix) <- colnames(across_cond_cor)
rownames(accuracy_matrix) <- colnames(across_cond_cor)

pheatmap::pheatmap(accuracy_matrix,
                   main ='binary classification accuracy >0.57 significant',
                   breaks =  seq(0, 1, by = 0.01))

# comparing differential covariation
signficant_atleast_once = DT_normalised[,.(max_bicor = max(bicor,na.rm = T)), by = .(Protein_1,Protein_2)
                                        ][max_bicor>0.57][,max_bicor:=NULL]
significant_pairs = merge(DT_normalised,
                          signficant_atleast_once, by = c('Protein_1','Protein_2'))


across_cond_cor = significant_pairs |> 
  dcast(Protein_1 + Protein_2 ~ Condition, value.var = 'condition_specific_cor') 
across_cond_cor |> dim()
across_cond_cor = f_BIC(across_cond_cor[,-c(1:2)],value_name = 'condition_specific_cor')

across_cond_cor_rev = copy(across_cond_cor)
across_cond_cor_rev[,`:=`(Protein_1=Protein_2,
                          Protein_2=Protein_1)]
across_cond_cor = rbind(across_cond_cor,
                        across_cond_cor_rev)
across_cond_cor_matrix = dcast(across_cond_cor,Protein_1 ~ Protein_2, 
                               value.var = 'condition_specific_cor') |> 
  tibble::column_to_rownames('Protein_1') |> 
  as.matrix()
diag(across_cond_cor_matrix) <- 1
pheatmap::pheatmap(across_cond_cor_matrix, 
                   main ='at least once significant condition specific covariation')




DT_normalised$zscore |> hist()
DT_normalised[zscore>4] |> View()

# DT_normalised[Protein_1 == 'A0A1D8PKE9' & 
#                 Protein_2 == 'A0A1D8PEK7'] |> View()
check_Protein_1 = 'Q5A0L0'
check_Protein_2 = 'A0A1D8PLA7'

top_covariation_parnters_1 = DT_normalised[Protein_1  == check_Protein_1| 
                                             Protein_2  == check_Protein_1
                                           ][,Partner:= fifelse(Protein_1 ==check_Protein_1,
                                                                Protein_2  ,Protein_1)]
top_covariation_parnters_1[,avg_bicor:=mean(bicor), by = .(Protein_1,Protein_2)]
top_covariation_parnters_1 |> ggplot(aes(x= avg_bicor,y = bicor))+
  geom_hex()+
  facet_wrap('Condition')

top_covariation_parnters_2 = DT_normalised[Protein_1  == check_Protein_2| 
                                             Protein_2  == check_Protein_2
][,Partner:= fifelse(Protein_1 ==check_Protein_2,
                     Protein_2  ,Protein_1)]
top_covariation_parnters_2[,avg_bicor:=mean(bicor), by = .(Protein_1,Protein_2)]
top_covariation_parnters_2 |> ggplot(aes(x= avg_bicor,y = bicor))+
  geom_hex()+
  facet_wrap('Condition')
top_partners_1 = top_covariation_parnters_1[order(-avg_bicor),.(Partner,avg_bicor)] |> unique()
top_partners_1 = top_partners_1[,Partner] |> head(15)
top_partners_2 = top_covariation_parnters_2[order(-avg_bicor),.(Partner,avg_bicor)] |> unique()
top_partners_2 = top_partners_2[,Partner] |> head(15)

df |> dim()

bicorcor_OI_int_piv  = df |>
  as.data.frame() |> 
  tibble::rownames_to_column('Uniprot') |> 
  as.data.table() |> 
  melt(id.vars = 'Uniprot', variable.name = 'experiment',value.name = 'LogRatio')

bicorcor_OI_int_piv[,condition:= str_extract(experiment,'_S[:print:]*') |> str_remove('^_')]
bicorcor_OI_int_piv[,category:=dplyr::case_when(
  Uniprot == check_Protein_1  ~check_Protein_1,
  Uniprot == check_Protein_2  ~check_Protein_2,
  Uniprot %in% top_partners_1 ~glue::glue('top_partners_{check_Protein_1}'),
  Uniprot %in% top_partners_2 ~glue::glue('top_partners_{check_Protein_2}'),
  TRUE~'other'
)]
bicorcor_OI_int_piv =bicorcor_OI_int_piv[category!='other']
bicorcor_OI_int_piv[,avg_prot_expression_condition := median(LogRatio,na.rm = T), by =  condition]
bicorcor_OI_int_piv[,condition_specific_expression:= LogRatio -avg_prot_expression_condition ]

conditions_to_plot = DT_normalised[Protein_1 == check_Protein_1 & 
                Protein_2 ==check_Protein_2][order(-bicor)]
conditions_to_plot = conditions_to_plot[c(1,nrow(conditions_to_plot)),
                                        Condition]

colour_pallet = c('#ffcccb','darkred','darkblue','lightblue')
names(colour_pallet) = c(glue::glue('top_partners_{check_Protein_1}'),check_Protein_1,
                         check_Protein_2,glue::glue('top_partners_{check_Protein_2}'))
bicorcor_OI_int_piv = bicorcor_OI_int_piv[category != 'other'
                                          ][condition %in% conditions_to_plot]
both_present = bicorcor_OI_int_piv[Uniprot %in% c(check_Protein_1,check_Protein_2) & 
                      !is.na(condition_specific_expression)
                    ][,.(both_prots_present = .N), by = experiment
                      ][both_prots_present==2,experiment]
bicorcor_OI_int_piv = bicorcor_OI_int_piv[experiment %in% both_present ]
  ggplot(bicorcor_OI_int_piv,aes(x = experiment, y= condition_specific_expression,
                                 colour= category,group=Uniprot))+
  # geom_line(data = bicorcor_OI_int_piv[category == 'other'],alpha = 0.05)+
  # geom_line()+
  geom_line(data = bicorcor_OI_int_piv[!(category %in% c(check_Protein_1,check_Protein_2))],alpha = 0.6)+
  geom_line(data = bicorcor_OI_int_piv[(category %in% c(check_Protein_1,check_Protein_2))],alpha = 1)+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_colour_manual(values = colour_pallet)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  facet_wrap('condition',nrow = 2)

pairwise_cor = df[rownames(df) %in% c(check_Protein_1,check_Protein_2),
                  str_detect(colnames(df),as.character(conditions_to_plot[1]))] |> t() |> 
  as.data.frame()
pairwise_cor |> ggplot(aes(x = .data[[check_Protein_1]] , y = .data[[check_Protein_2]]  ))+
  geom_point()+
  geom_abline(intercept = 0, slope =1)+
  ggtitle('bicor doesnt capture similar behavior in partitions',
          subtitle = 'when the median normalisation happens globally')

# calculating_groups_distances <- function(){
#   partners = bicorcor_OI_int_piv[str_detect(category,'top_partners')] |> 
#     dcast(Uniprot~)
  
  DT_partners = DT_normalised[(Protein_1 %in% c(check_Protein_1,
                                                check_Protein_2,
                                                top_partners_1,
                                                top_partners_2) &
                  Protein_2 %in% c(check_Protein_1,
                                   check_Protein_2,
                                   top_partners_1,
                                   top_partners_2)) 
                & Condition == as.character(conditions_to_plot[2])]  
  DT_partners_rev = copy(DT_partners)
  DT_partners_rev[,`:=`(Protein_1=Protein_2,
                      Protein_2=Protein_1)]
  DT_partners = rbind(DT_partners,
                      DT_partners_rev)
  DT_partners_matrix = dcast(DT_partners,Protein_1 ~ Protein_2, value.var = 'bicor') |> 
    tibble::column_to_rownames('Protein_1') |> 
    as.matrix()
  diag(DT_partners_matrix) <- 1
  annotation_col = bicorcor_OI_int_piv[,.(Uniprot,category)] |> 
    unique() |> 
    tibble::column_to_rownames('Uniprot')
  DT_partners_matrix[is.na(DT_partners_matrix)] <- 0
  pheatmap::pheatmap(DT_partners_matrix,annotation_col = annotation_col,
                     annotation_row = annotation_col, main = conditions_to_plot[2] )
# }

DT_normalised[Protein_1 == check_Protein_1 & Protein_2 %in% top_partners_2|
                Protein_2 == check_Protein_1 & Protein_1 %in% top_partners_2] |> 
  ggplot(aes(x = Condition, y = zscore))+
  geom_boxplot()


DT_avg = DT_normalised[!is.na(bicor),.(N_covariations = .N, 
      avg_covariation = mean(bicor, na.rm = T),
      sd_covariation = sd(bicor, na.rm = T)), by = .(Protein_1,Protein_2)]

DT_avg[avg_covariation>0.3][order(-avg_covariation)]|> 
  ggplot(aes(x= avg_covariation, y = sd_covariation))+
  geom_hex()+
  scico::scale_fill_scico()+
  ggtitle('High covarying pairs, covary similarly across conditions')
ggsave(here::here('out','plots','hexplot_covaration_avg_conditions.pdf'))

DT_normalised[,.SD[sample(.N,1000000)], by = .(Condition)] |> ggplot(aes(x = avg_cor, y = bicor))+
  geom_hex()+
  theme_bw()+
  facet_wrap('Condition')+
  annotate("rect", xmin = 0.57, xmax = 1, ymin = -1, ymax = 1,
           alpha = .2, colour = 'darkred')+
  annotate("rect", xmin = -1, xmax = 1, ymin = 0.57, ymax = 1,
           alpha = .2, colour = 'darkred')
ggsave(here::here('out','plots','hexplot_conditions_similarity_sample.pdf'))
