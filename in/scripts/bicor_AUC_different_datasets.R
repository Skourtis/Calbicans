# check bicor against GS
library(data.table);library(ggplot2);library(WGCNA)
library(PRROC);library(ggpubr);library(distanceHD)
# Uniprot mapping
Uniprot_annot = fread(here::here('in','datasets','uniprotkb_taxonomy_id_237561_2025_07_08.tsv.gz'))
Uniprot_annot[,CGD:=str_remove(CGD,';[:print:]*$')]

# GoldStandard
GS = fread(here::here('out','datasets','STRING_TP_FP.gz') )
#### Annotate with gold standard ####
setnames(GS,c('Protein_1','Protein_2','Class'))



# load datasets
KO= fread(here::here('in','datasets','kinase_ko_proteomes.tsv'))
KO = KO[str_detect(Strain,'^CAL')
        ][!is.na(PG.MaxLFQ)][,PG.MaxLFQ:=log2(PG.MaxLFQ)]
KO[,mean_expression:=median(PG.MaxLFQ,na.rm = T),by = .(Strain,replicate)]
KO[,PG.MaxLFQ := PG.MaxLFQ -mean_expression ]
  ggplot(KO,aes(x = Strain,colour = as.factor(replicate),y = PG.MaxLFQ))+
  geom_boxplot()
KO[,Strain:=paste0('strain_',Strain)]

KO_avg =  KO[,.(abundance = mean(PG.MaxLFQ)),by = .(Genes,Strain)]

  KO_avg =  KO_avg |> dcast(Genes~Strain,value.var = 'abundance') 
  KO_avg |> nrow()
  KO_avg = merge(KO_avg,Uniprot_annot[,.(Entry,CGD)],by.x = 'Genes',by.y = 'CGD')
  KO_avg |> nrow() 
  KO_avg[,Genes:=NULL]
  KO_avg <- KO_avg |> tibble::column_to_rownames('Entry')
  KO_avg |> dim()
  
  KO_rep1 =  KO[replicate==1] |> dcast(Genes~Strain,value.var = 'PG.MaxLFQ') 
  KO_rep1 |> nrow()
  KO_rep1 = merge(KO_rep1,Uniprot_annot[,.(Entry,CGD)],by.x = 'Genes',by.y = 'CGD')
  KO_rep1 |> nrow() 
  KO_rep1[,Genes:=NULL]
  KO_rep1 <- KO_rep1 |> tibble::column_to_rownames('Entry')

  # 
  KO_rep2 =  KO[replicate==2] |> dcast(Genes~Strain,value.var = 'PG.MaxLFQ') 
  KO_rep2 |> nrow()
  KO_rep2 = merge(KO_rep2,Uniprot_annot[,.(Entry,CGD)],by.x = 'Genes',by.y = 'CGD')
  KO_rep2 |> nrow() 
  KO_rep2[,Genes:=NULL]
  KO_rep2 <- KO_rep2 |> tibble::column_to_rownames('Entry')
  
  # wild_type cells
  
  WT= fread(here::here('in','datasets','updated_synergy_proteomes.tsv'))
  WT = WT[!is.na(CGDID)]
  WT[,.(.N),by = .(CGDID,plate_position)][,N] |> hist()
  WT[,abundance:= log2(abundance)]
  WT[,mean_expression:=median(abundance,na.rm = T),by = .(plate_position)]
  WT[,abundance := abundance -mean_expression ]

  WT = WT|> dcast(CGDID~plate_position,value.var = 'abundance') 
  WT |> nrow()
  WT = merge(WT,Uniprot_annot[,.(Entry,CGD)],by.x = 'CGDID',by.y = 'CGD')
  WT |> nrow() 
  WT[,CGDID:=NULL]
  WT <- WT |> tibble::column_to_rownames('Entry')
  WT |> dim()
  
  # creating a subset to match the size of the KO dataset
  WT_83 <- WT[,sample(1:ncol(WT),83)]
  
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

  subset_stress_83 <- subset_stress[,sample(1:ncol(subset_stress),83)]
  subset_stress_83 |> dim()
  
  list_dataset <- list(KO_avg = KO_avg,
                       KO_rep1 = KO_rep1,
                       KO_rep2 = KO_rep2,
                       WT_isolates = WT,
                       WT_isolates_83 = WT_83,
                       stressed_isolates = subset_stress,
                       subset_stress_83 = subset_stress_83)
  
  
  
list_PRcurves = list()
n_detected = 30
for(dataset in names(list_dataset)){

df = list_dataset[[dataset]]
# Select only proteins with at least ratios
feature_count <- apply(df, 1, function(x){ sum(!is.na(x))} )
# almost all proteins detected in all samples
feature_count |> hist()
df <- df[feature_count >= n_detected,]

# Transpose and median normalise the rows to avoid spurious correlations
df <- t(df)
tmp_medians <- apply( df, 1, median, na.rm = TRUE )  
df <- sweep( df, 1, tmp_medians, FUN = "-" )


#### Calculate correlations between proteins ####

# Define function for bicor correlation. It takes an input data frame and a name of the output value. 
# It outputs a row-wise sorted list of protein pairs along with calculated distance metric
f_BIC <- function( input_data, value_name ){
  missingness_tmp <- base::crossprod(!is.na(input_data))>10 # with microproteins being very rare it's common to only have 1 point in common produce fake correlations
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
list_PRcurves[[dataset]] <- f_BIC(df, dataset)

}
# Combine the results
DT <- Reduce(function(...) merge(...,by = c("Protein_1", "Protein_2")),list_PRcurves)    # This works because all pairs are sorted A > B
DT[,avg_covariation:=mean(c(KO_rep1,KO_rep2),na.rm = T),by = .(Protein_1,Protein_2)]
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
DT[, random_classifier := sample(`KO_avg`)]


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

# Set standard plotting parameters
plot_cols <- c( 
  # prothd2_prohd1_subset = "#37abc8ff", 
                KO = "#ff2a7fff", random_classifier = "grey50") # Set plotting colours
plot_labels <- c(
  # prothd2_prohd1_subset = "proHD2_subset", 
                 KO = "KO", random_classifier = "Random")

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
ggsave(here::here('out','plots','STRING_GS_bicor_PRcurves.pdf'),p)

# make a plot like proHD2 with the different sizes, and the total proteomes they are able to characterise
# check size and MAD per protein per dataset

dataset_size = data.table(Datasets = c('WT isolates','Stressed Isolates','Kinase KO'),
                        N_Prot = c(nrow(WT),nrow(subset_stress),nrow(KO_avg)),
                        conditions = c(ncol(WT),ncol(subset_stress),ncol(KO_avg)))
dataset_sizes = ggplot() + 
  geom_rect(aes(xmin = 0, xmax = unlist(dataset_size[2,3]), ymin = 0, ymax = unlist(dataset_size[2,2])),fill = "#7E2811", colour = 'black', alpha = 0.8)+
  geom_rect(aes(xmin = 0, xmax = unlist(dataset_size[1,3]), ymin = 0, ymax = unlist(dataset_size[1,2])), fill = "grey50", colour = 'black',alpha = 0.8)+
  geom_rect(aes(xmin = 0, xmax = unlist(dataset_size[3,3]), ymin = 0, ymax = unlist(dataset_size[3,2])), fill = "darkblue", colour = 'black',alpha = 0.8)+
  annotate("text", x = unlist(dataset_size[1,3])/2, y = unlist(dataset_size[1,2])+300, label = dataset_size[1,Datasets],  size= 3)+
  annotate("text", x = unlist(dataset_size[2,3])/2, y = unlist(dataset_size[2,2])+300, label = dataset_size[2,Datasets],  size= 3)+
  annotate("text", x = unlist(dataset_size[3,3])/2, y = unlist(dataset_size[3,2])+300, label = dataset_size[3,Datasets],  size= 3)+
  theme_bw()+
  xlab('number of conditions')+
  ylab('Number of proteins')
ggsave(here::here('out','plots','dataset_sizes.pdf'),dataset_sizes)

# as well as MAD per TP and FP 

common_proteins = Reduce(intersect,list(rownames(WT),
                       rownames(subset_stress),
                       rownames(KO_avg))
                       )
WT_mad = matrixStats::rowMads(as.matrix(WT),na.rm = T) |> 
  tibble::enframe(value = 'mad',name = 'Uniprot') |> 
  as.data.table()

KO_mad = matrixStats::rowMads(as.matrix(KO_avg),na.rm = T) |> 
  tibble::enframe(value = 'mad',name = 'Uniprot') |> 
  as.data.table()

subset_stress_mad = matrixStats::rowMads(as.matrix(subset_stress),na.rm = T) |> 
  tibble::enframe(value = 'mad',name = 'Uniprot') |> 
  as.data.table()

mad_datasets = Reduce(rbind,list(WT_mad[,dataset:='WT_isolates'],
                                 subset_stress_mad[,dataset:='stressed_isolates'],
                                 KO_mad[,dataset:='kinase_KO']
))
mad_datasets[,Common_proteins:= fifelse(Uniprot %in% common_proteins, 'common','unique')]

stat_box_data <- function(y) {
  return( 
    data.frame(
      y = 0.5+1.1*max(y),  # may need to modify this depending on your data
      label = paste('n =', length(y))
    )
  )
}




mad_datasets |> ggplot(aes(x =dataset,y = mad , fill =  dataset))+
  geom_boxplot()+
  theme_bw()+
  facet_wrap('Common_proteins')+
  stat_summary(
    fun.data = stat_box_data, 
    geom = "text", 
    hjust = 0.5,
    vjust = 0.9
  )+
  theme(legend.position = 'bottom')
ggsave(here::here('out','plots','dataset_mad.pdf'))

WT_missingness = matrixStats::rowMeans2(as.matrix(is.na(WT))) |> 
  tibble::enframe(value = 'Missingness',name = 'Uniprot') |> 
  as.data.table()

KO_missingness = matrixStats::rowMeans2(as.matrix(is.na(KO_avg))) |> 
  tibble::enframe(value = 'Missingness',name = 'Uniprot') |> 
  as.data.table()

subset_stress_missingness = matrixStats::rowMeans2(as.matrix(is.na(subset_stress))) |> 
  tibble::enframe(value = 'Missingness',name = 'Uniprot') |> 
  as.data.table()

missingness_datasets = Reduce(rbind,list(WT_missingness[,dataset:='WT_isolates'],
                                         subset_stress_missingness[,dataset:='stressed_isolates'],
                                         KO_missingness[,dataset:='kinase_KO']
))
missingness_datasets[,Common_proteins:= fifelse(Uniprot %in% common_proteins, 'common','unique')]

missingness_datasets |> ggplot(aes(x =dataset,y = Missingness , fill =  dataset))+
  geom_boxplot()+
  theme_bw()+
  facet_wrap('Common_proteins')+
  stat_summary(
    fun.data = stat_box_data, 
    geom = "text", 
    hjust = 0.5,
    vjust = 0.9
  )+
  theme(legend.position = 'bottom')
ggsave(here::here('out','plots','dataset_missingness.pdf'))


GS_mad = merge(GS,subset_stress_mad, by.x = 'Protein_1', by.y = 'Uniprot' ) |> 
  merge(subset_stress_mad, by.x = 'Protein_2', by.y = 'Uniprot' )

GS_mad[,.SD[sample(.N,10000)],by = Class] |> 
  ggplot(aes (x= mad.x, y =mad.y, colour = Class ))+
  geom_density2d(linemitre = 1)+
  labs(x = 'Protein 1 MAD',y = 'Protein 2 MAD')+
  # theme_bw()+
  ggtitle('mad of protein pairs TP/FP')+
  theme_bw()

# ggsave(here::here('out','plots','mad_GS.pdf'))
# +facet_wrap('GS')


seq(min(subset_stress_mad$mad),
    max(subset_stress_mad$mad),length.out = 10)
dplyr::ntile(subset_stress_mad$mad,3) |> table()
GS_mad[,`:=`(mad1 = dplyr::ntile(mad.x,5),
               mad2 = dplyr::ntile(mad.y,5))]
GS_mad[,bands:= paste(mad1,mad2,sep = '_')]
set.seed(123)
GS_mad[Class=='TP',N_TP := .N, by =bands ]
total_TP = nrow(GS_mad[Class=='TP'])
GS_mad[,N_TP := mean(N_TP,na.rm = T), by =bands ]

GS_mad[Class=='TP',N_TP := .N, by =bands ]
GS_mad[,Perc_TP:= N_TP/total_TP]
total_FP = total_TP*15
sampled_GS = GS_mad[!is.na(N_TP),.SD[sample(.N,min(.N,total_FP*Perc_TP))],by = .(bands,Class)] 
sampled_GS = sampled_GS[,.(Protein_1,Protein_2,Class)]
sampled_GS[,id:=1:.N]
set.seed(1234)
train_set = sample(sampled_GS$id,nrow(sampled_GS)*0.7)
sampled_GS[,set:= fifelse(id %in% train_set,'train','test')]
fwrite(sampled_GS[,id:=NULL],here::here('out','datasets','sampled_GS.gz'))
new_GS = rbind(GS_mad[,dataset:='unfiltered'],
               sampled_GS[,dataset:='sampled'] )
new_GS[,.N,by = .(Class,dataset)]
new_GS[,.N,by = .(Class,dataset,bands)] |> 
  ggplot(aes(x = bands, y = log2(N), fill = Class))+
  geom_col(position = 'dodge')+
  facet_wrap('dataset',scale= 'free', nrow =2)
new_GS= new_GS[,.SD[sample(.N,total_TP)],by = .(Class,dataset)]
ggplot(new_GS,aes (x= mad.x, y =mad.y, colour = Class ))+
  geom_density_2d()+
  labs(x = 'Protein 1 MAD',y = 'Protein 2 MAD')+
  # theme_bw()+
  ggtitle('mad of protein pairs TP/FP')+
  theme_bw()+
  facet_wrap('dataset')
ggsave(here::here('out','plots','mad_GS.pdf'))
