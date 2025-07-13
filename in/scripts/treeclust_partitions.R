# stressed isolates does the best in bicor
# so we try to optimize a covariation strategy for it
library(data.table);library(ggplot2);library(treeClust)
library(PRROC);library(ggpubr);library(WGCNA);library(furrr)

# Uniprot mapping
Uniprot_annot = fread(here::here('in','datasets','uniprotkb_taxonomy_id_237561_2025_07_07.tsv.gz'))
Uniprot_annot[,CGD:=str_remove(CGD,';[:print:]*$')]

# Load the protein ratios 
setDTthreads(90)
# gold_standards
GS = fread(here::here('out','datasets','sampled_GS.gz'))
# Load the protein ratios 
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
subset_stress= subset_stress[!is.na(Genes)]
subset_stress |> nrow() 
subset_stress[,Genes:=NULL]
subset_stress <- subset_stress |> tibble::column_to_rownames('Entry')
subset_stress |> dim()

feature_count <- apply(subset_stress, 1, function(x){ sum(!is.na(x))} )
# remove entries with less than 40 datapoints
df = subset_stress[feature_count >=40,]
tmp_medians <- apply( df, 1, median, na.rm = TRUE )  
df <- sweep( df, 1, tmp_medians, FUN = "-" )
# Define function for bicor correlation. It takes an input data frame and a name of the output value. 
# It outputs a row-wise sorted list of protein pairs along with calculated distance metric
f_BIC<- function( input_data, value_name ){
  missingness_tmp <- base::crossprod(!is.na(input_data))>40 # with microproteins being very rare it's common to only have 1 point in common produce fake correlations
  #using suggested bicor
  tmp <- bicor( input_data , use = "pairwise.complete.obs", nThreads = 6)                     # Robust correlation using all pairwise complete observations
  tmp[missingness_tmp == FALSE] <- NaN
  tmp <- as.data.table( reshape2::melt( as.matrix( tmp  )))                                   # Turn distance matrix into a pair-wise data.table
  tmp <- tmp[, .( Protein_1 = as.character(Var1), Protein_2 = as.character(Var2), value ) ]   # Rename and change to character
  tmp <- tmp[ Protein_1 > Protein_2 ]                                                         # Remove redundant pairs by keeping only A > B, removing A == B and B < A pairs
  names(tmp)[3] <- value_name                                                                 # Assign new value name
  return(tmp)
}
f_PCC<- function( input_data, value_name ){
  missingness_tmp <- base::crossprod(!is.na(input_data))>10 # with microproteins being very rare it's common to only have 1 point in common produce fake correlations
  #using suggested bicor
  tmp <- cor( input_data , use = "pairwise.complete.obs", nThreads = 6)                     # Robust correlation using all pairwise complete observations
  tmp[missingness_tmp == FALSE] <- NaN
  tmp <- as.data.table( reshape2::melt( as.matrix( tmp  )))                                   # Turn distance matrix into a pair-wise data.table
  tmp <- tmp[, .( Protein_1 = as.character(Var1), Protein_2 = as.character(Var2), value ) ]   # Rename and change to character
  tmp <- tmp[ Protein_1 > Protein_2 ]                                                         # Remove redundant pairs by keeping only A > B, removing A == B and B < A pairs
  names(tmp)[3] <- value_name                                                                 # Assign new value name
  return(tmp)
}

# original trC parameter on whole ProHD2
f_trC <- function( input_data, cond_vector ){
  surrogatestyle_input = cond_vector[1]
  minsplits = cond_vector[2]
  cps = cond_vector[3] 
  serule = cond_vector[4] 
  partition= cond_vector[5]
  value_name = glue::glue('surrogate_{surrogatestyle_input}_splits_{minsplits}_cp_{cps}_serule_{serule}_partition_{partition}')
  
  # missingness_tmp <- base::crossprod(!is.na(input_data))>10 # with microproteins being very rare it's common to only have 1 point in common produce fake correlation
  input_modified_for_treeClust <- as.data.frame( t(input_data )   )                                      # Transpose and turn to data frame, as required by this function. 
  colnames(input_modified_for_treeClust) <- paste("feature", 1:ncol(input_modified_for_treeClust), sep = "_")     # ... also treeClust can't handle feature names starting with numbers apparently
  tmp = list()
  seed = 1234
  set.seed(seed)
  
  while(!is.data.table(tmp) & seed < 1250){
    seed = seed+1
    tryCatch({
      tmp <- treeClust.dist( input_modified_for_treeClust , d.num = 2,  verbose = T,
                             rcontrol = rpart.control(cp = cps, surrogatestyle  = surrogatestyle_input,
                                                      minsplit = minsplits),
                             control = treeClust.control(serule =  serule)
      )
      tmp = as.matrix(tmp)
      # tmp[missingness_tmp == FALSE] <- NaN
      tmp <- as.data.table( reshape2::melt( as.matrix( tmp  )))
      tmp <- tmp[, .( Protein_1 = as.character(Var1), Protein_2 = as.character(Var2), value = (1-value)) ]
      setkeyv(tmp,c('Protein_1','Protein_2'))
      tmp <- tmp[ Protein_1 > Protein_2 ]
      names(tmp)[3] <- glue::glue('seed_{seed}_{value_name}')
      
    },
    error=function(e) {
      message('An Error Occurred')
      print(e)
    },
    warning=function(w) {
      message('A Warning Occurred')
      print(w)
    }
    )}
  if(is.data.table(tmp)){
    return(tmp)
  }else{
    print(seed)
    return('error')
  }
}
bicor_corr = f_BIC(t(df),'corr_bicor')
pcc_corr = f_PCC(t(df),'corr_PCC')

# original ProHD parameters
tree_clust_corr = f_trC(t(df),c(surrogatestyle_input = 0,
                                minsplits = 20,
                                cps = 0.01,
                                serule =0,
                                partition = 20))
Whole_ProHD2 = Reduce(function(...) merge(...,by = c("Protein_1", "Protein_2")),list(ProHD2_bicor = bicor_corr,
                                                                                     ProHD2_PCC = pcc_corr,
                                                                                     ProHD2_tree = tree_clust_corr))
fwrite(Whole_ProHD2, here::here('out','datasets','whole_stressed_covariations.gz'))    
fwrite(bicor_corr, here::here('out','datasets','bicor_all_samples.gz'))    

Whole_ProHD2 = fread(here::here('out','datasets','whole_stressed_covariations.gz'))


df |> dim()
df[1:2,1:2]
cut(colnames(df),5)
partition_sizes= list(partition_1 = 1:200,
                 partition_2 = 201:400,
                 partition_3 = 401:600,
                 partition_4 = 601:800,
                 partition_5 = 801:1030)
partition_list = purrr::map(.x =partition_sizes, ~t(df[,.x]) )

set.seed(1234)
GS_train = GS[set =='train']
plan(sequential)
plan(multisession, workers = 20)

# grid searching AUC for each partition
# looping through each of the parameters
complete_AUCs = data.table()
for(partition in 1:length(partition_list)){
  # define grid
  surrogates  = c(0,1)
  minsplits = 20
  cps = 0.015
  serule = c(0,0.5,1,1.5,1.8,2,2.5,3.5,4,8)
  
  conditions = expand.grid(surrogates,minsplits,cps, serule)
  count = 0
  while(count <3){
    print(count)
    if(count == 0){
      name_hyper = 'surrogate_serule'
    }
    if(count == 1){
      conditions = expand.grid(surrogates,minsplits,cps, serule)
      minsplits = c(10,20,50,100,300,400,500,600,1000)
      conditions = expand.grid(surrogates,minsplits,cps, serule)
      name_hyper = 'surrogate_serule_minplit'
    }
    if(count == 2){
      
      cps= c(0.1,0.01,0.0125,0.015,0.001,0.0075,0.0001,0.00001)
      conditions = expand.grid(surrogates,minsplits,cps, serule)
      name_hyper = 'surrogate_serule_minplit_cps'
    }
    print(partition)
    partition = as.numeric(partition)
    conditions$partition = as.numeric(partition)
    conditions_to_use = apply(conditions,1,as.vector, simplify = F)
    all_partition_list = rep(partition_list[partition],each =length(conditions_to_use)) 
    partition= as.character(partition)
    DT_partition = data.table()
    
    # hyper_param_list_ = purrr::map2(.x = all_partition_list[1],.y = conditions_to_use[1],  ~f_trC(.x, .y))
    
    hyper_param_list_ = furrr::future_map2(.x = all_partition_list,.y = conditions_to_use,  ~f_trC(.x, .y))
    hyper_param_list_ = purrr::keep(hyper_param_list_,is.data.table)
    DT_partition = Reduce(function(...) merge(...,by = c("Protein_1", "Protein_2")),hyper_param_list_)
    rm(hyper_param_list_)
    fwrite(DT_partition,here::here('out','datasets',glue::glue('hyperparam_partition__bicor_{partition}_{name_hyper}_balanced.gz')))
    
    DT_partition <- DT_partition[ Protein_1 != Protein_2 ]                  # Removing self-associations that arose between isoforms (can't test those anyway)
    DT_partition[ Protein_1 <= Protein_2 ]                        # Validate that A > B still fullfilled   
    standard = 'train'
    GS_train = GS[set == standard][,set := NULL ]
    setnames(GS_train,c('Protein_1','Protein_2','Class'))
    
    GS_train[ Protein_1 <= Protein_2, ]     # All pairs sorted A > B, so ready to merge with DT 
    DT <- merge(DT_partition, GS_train, by = c("Protein_1", "Protein_2"))
    DT[, .N, Class]  # Currently about 3% TP pairs
    
    # Downsample FPs to 90% - this is essential to make curves comparable across project
    # n_FP_to_remove <- DT[ Class == "FP", .N ] - DT[ Class == "TP", .N ] * 9
    # DT[ sample( which( Class == "FP" ), n_FP_to_remove ), Class := NA ]    # From the rows where the class is FP, randomly sample a certain subset and turn to NA
    DT <- DT[ !is.na(Class)  ][Class !='']                                             # Remove those rows
    DT[, .N, Class]  # Now 10% TP pairs
    
    # Add a randomised classifier to the input table
    sample_from = DT[,3] |> unlist()
    DT[, random_classifier := sample(sample_from)]
    
    
    #### Run precision recall analyses ####
    
    #### Define PR curve function and PR input data sets ####
    
    # An R function for precision - recall analysis using the PRROC package.
    # Inputs: data = a data table to analyse; gs = the column containing the class labels; metric = the column containing the metric to be tested
    # Optional input: the "modifier" argument specifies how pairs will be ranked. The default "descending" means that pairs are sorted such that higher values correspond to higher confidence of association. For other metrics (e.g. distances), set to "ascending" to reverse the order, or to "absolute" to rank pairs by the absolute metric value
    # Output: A list of length 2, containing some stats and the actual curve coordinates, respectively
    PR_curve <- function(data, gs, metric, modifier = "descending"){
      # data = DT
      # gs = 'Class'
      # metric = name
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
      tmp_curve_dt <- rbind( tmp_curve_dt[ Recall  < 0.05                ][ sample(.N, min(.N,180)) ],   # Downsample the curves for faster plotting. Keep more points for the low recall areas because they are less robust.
                             tmp_curve_dt[ Recall >= 0.05 & Recall < 0.1 ][ sample(.N,  min(.N,200)) ],
                             tmp_curve_dt[ Recall >= 0.10 & Recall < 0.2 ][ sample(.N, min(.N,200)) ],
                             tmp_curve_dt[ Recall >= 0.20 & Recall < 0.4 ][ sample(.N, min(.N,200)) ],
                             tmp_curve_dt[ Recall >= 0.40                ][ sample(.N, min(.N,200)) ])
      output <- list( 
        "stats" = data.table( "AUPRC" = tmp_curve$auc.integral, "AUPRC_random" = tmp_curve$rand$auc.integral, "N_TP" = length(tmp_TP), gold_standard = gs, metric = metric), 
        "curve" = tmp_curve_dt)
      return(output)
    }
    
    # Perform the PR analyses using my custom function. 
    PR_curves = list()
    PR_stats = data.table()
    PR_data =data.table()
    for(name in colnames(DT)[-c(1:2,ncol(DT)-1)]){
      print(name)
      PR_curves[[name]] = PR_curve( DT, "Class", name)
      PR_stats = rbind(PR_stats,
                       PR_curves[[name]][['stats']])
      PR_data = rbind(PR_data,
                      PR_curves[[name]][['curve']])
    }
    PR_stats[,partition:=partition]
    PR_stats[,set:=standard]
    PR_stats= PR_stats[metric!= 'random_classifier']
    # AUCs = rbind(AUCs, PR_stats)
    complete_AUCs = rbind(complete_AUCs, PR_stats)
    fwrite(complete_AUCs,here::here('out','datasets', 'complete_AUCs_balanced.csv'))
    PR_stats = PR_stats[order(-AUPRC)] |> head(7)
    if(count == 0){
      surrogates = PR_stats$metric |> stringr::str_match('surrogate_(.)_')
      surrogates = surrogates[,2] |> unique() |> as.numeric()
      serule = PR_stats$metric |> stringr::str_match('serule_([:print:]*?)_')
      serule = serule[,2] |> unique() |> as.numeric()
    }
    if(count == 1){
      surrogates = PR_stats$metric |> stringr::str_match('surrogate_(.)_')
      surrogates = surrogates[,2] |> unique() |> as.numeric()
      serule = PR_stats$metric |> stringr::str_match('serule_([:print:]*?)_')
      serule = serule[,2] |> unique() |> as.numeric()
      minsplits = PR_stats$metric |> stringr::str_match('splits_([:print:]*?)_')
      minsplits = minsplits[,2] |> unique() |> as.numeric()
    }
    if(count == 2){
      surrogates = PR_stats$metric |> stringr::str_match('surrogate_(.)_')
      surrogates = surrogates[,2] |> unique() |> as.numeric()
      serule = PR_stats$metric |> stringr::str_match('serule_([:print:]*?)_')
      serule = serule[,2] |> unique() |> as.numeric()
      minsplits = PR_stats$metric |> stringr::str_match('splits_([:print:]*?)_')
      minsplits = minsplits[,2] |> unique() |> as.numeric()
      cps = PR_stats$metric |> stringr::str_match('cp_([:print:]*?)_')
      cps = cps[,2] |> unique() |> as.numeric()
    }
    count = count+1
  }
}
# complete AUCs for all partitions
complete_AUCs = fread(here::here('out','datasets', 'complete_AUCs_balanced.csv'))
min_scale= complete_AUCs$AUPRC_random |> min()
max_scale = complete_AUCs$AUPRC |> max()
complete_AUCs |> 
  ggplot(aes(x = AUPRC_random, y = AUPRC, colour = partition))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  lims(x= c(min_scale,max_scale),
       y= c(min_scale,max_scale))
# best hyperparame is much better than random
hyperparams = complete_AUCs[order(-AUPRC)][,head(.SD,1), by = partition][,.(partition,metric)]
hyper_paramlist = list()
bicor_partition = data.table()
for(i in 1:length(partition_list)){
  metric_tmp = hyperparams[partition ==i,metric]
  hyper_paramlist[[i]] =fread(here::here(glue::glue(here::here('out','datasets','hyperparam_partition__bicor_{i}_surrogate_serule_minplit_cps_balanced.gz'))),
                              select = c('Protein_1','Protein_2',metric_tmp))
  # bicor_partition = rbind(bicor_partition,f_BIC(partition_list[[i]], 'bicor_partition'))
  
}

# checkign how to combine the different partition scores
quality_partition = complete_AUCs[order(-AUPRC)][,head(.SD,1), by = partition][,Quality := AUPRC - AUPRC_random]
quality_partition = quality_partition[,.(partition)][,rank:= 1:.N]
hyper_paramlist_2 = hyper_paramlist[1:length(partition_list)]
hyper_paramlist_2 = purrr::map(.x = hyper_paramlist_2,~dplyr::rename(.x,'corr' = colnames(.x)[3]))
# best_AUCS_list = best_AUCS_list_comb[!is.na(corr)][,.(N_points=.N,
#                                                       min_treeclust_cor = min(corr,na.rm = T),
#                                                       mean_cor = mean(corr,na.rm = T),
#                                                       median_cor = median(corr,na.rm = T)
# ),     by = .(Protein_1,Protein_2)]

scaled_df = Reduce(function(...) merge(...,all = T,by = c("Protein_1", "Protein_2")),hyper_paramlist[1:length(partition_list)])
scaled_df_quantile =cbind(scaled_df[,1:2],
                          scaled_df[,3:7] |>  broman::normalize()) 
# IQR the scores
fwrite(scaled_df_quantile |> as.data.table(),here::here('out','datasets', 'partitions_IQR_scale_balanced.gz'))
scaled_df_quantile = fread(here::here('out','datasets', 'partitions_IQR_scale_balanced.gz'))
scale_values <- function(x){(x-min(x,na.rm = T))/(max(x,na.rm=T)-min(x,na.rm = T))}

scaled_df_top3 = scaled_df_quantile |> melt(id.vars = c('Protein_1','Protein_2'))
scaled_df_top3 = scaled_df_top3[!is.na(value)][order(-value)]
setkeyv(scaled_df_top3,c('Protein_1','Protein_2'))
scaled_df_top3 = scaled_df_top3[,.SD[head(.N,min(.N,3))], by = .(Protein_1,Protein_2)]
scaled_df_top3 = scaled_df_top3[,.(top_3_scaled= mean(value, na.rm = T)), by = .(Protein_1,Protein_2)]
# scaled_df_top3 = cbind(scaled_df_quantile[,1:2],
fwrite(scaled_df_top3, here::here('out','datasets', 'scaled_top_3_balanced.gz'))
# best_AUCS_list[,N_points:=NULL]
# hyper_paramlist[[11]] = best_AUCS_list

# compare test and train dataset
PR_stats_set = data.table()
hyper_paramlist_2 = Reduce(function(...) merge(...,by = c("Protein_1", "Protein_2")),list(
  scaled_df_quantile,
  scaled_df_top3
  ))
hyper_paramlist_2 <- hyper_paramlist_2[ Protein_1 != Protein_2 ]                  # Removing self-associations that arose between isoforms (can't test those anyway)
hyper_paramlist_2[ Protein_1 <= Protein_2 ]                        # Validate that A > B still fullfilled   
# best_AUCS_list[,.(Protein_1,Protein_2)]

plots_list = list()
for(standard in c('test','train')){
  print(standard)
  
  GS_tmp = GS[set == standard][,set := NULL ]
  setnames(GS_tmp,c('Protein_1','Protein_2','Class'))
  
  GS_tmp[ Protein_1 <= Protein_2, ]     # All pairs sorted A > B, so ready to merge with DT 
  # # 
  # # GS = GS_all_original[, c('OLN_1','OLN_2',gold_standard), with = F]
  # # setnames(GS,c('Protein_1','Protein_2','Class'))
  # # 
  # GS[ Protein_1 <= Protein_2, ]     # All pairs sorted A > B, so ready to merge with DT 
  
  # Annotate correlations with gold standards by merging
  DT <- merge(hyper_paramlist_2, GS_tmp, by = c("Protein_1", "Protein_2"))
  DT[, .N, Class]  # Currently about 3% TP pairs
  
  # Downsample FPs to 90% - this is essential to make curves comparable across project
  # n_FP_to_remove <- DT[ Class == "FP", .N ] - DT[ Class == "TP", .N ] * 9
  # DT[ sample( which( Class == "FP" ), n_FP_to_remove ), Class := NA ]    # From the rows where the class is FP, randomly sample a certain subset and turn to NA
  DT <- DT[ !is.na(Class)  ][Class !='']                                             # Remove those rows
  DT[, .N, Class]  # Now 10% TP pairs
  
  # Add a randomised classifier to the input table
  sample_from = DT[,3] |> unlist()
  DT[, random_classifier := sample(sample_from)]
  
  
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
    tmp_curve_dt <- rbind( tmp_curve_dt[ Recall  < 0.05                ][sample(min(.N, 180)) ],   # Downsample the curves for faster plotting. Keep more points for the low recall areas because they are less robust.
                           tmp_curve_dt[ Recall >= 0.05 & Recall < 0.1 ][ sample(min(.N, 200)) ],
                           tmp_curve_dt[ Recall >= 0.10 & Recall < 0.2 ][ sample(min(.N, 200))  ],
                           tmp_curve_dt[ Recall >= 0.20 & Recall < 0.4 ][ sample(min(.N, 200))  ],
                           tmp_curve_dt[ Recall >= 0.40                ][ sample(min(.N, 200)) ])
    output <- list( 
      "stats" = data.table( "AUPRC" = tmp_curve$auc.integral, "AUPRC_random" = tmp_curve$rand$auc.integral, "N_TP" = length(tmp_TP), gold_standard = gs, metric = metric), 
      "curve" = tmp_curve_dt)
    return(output)
  }
  
  # Perform the PR analyses using my custom function. 
  PR_curves = list()
  PR_stats = data.table()
  PR_data =data.table()
  for(name in colnames(DT)[-c(1:2,ncol(DT)-1)]){
    print(name)
    PR_curves[[name]] = PR_curve( DT, "Class", name)
    PR_stats = rbind(PR_stats,
                     PR_curves[[name]][['stats']])
    PR_data = rbind(PR_data,
                    PR_curves[[name]][['curve']])
  }
  PR_stats_set  = rbind(PR_stats_set,PR_stats[,set:= standard])
  
  #### Plot the PR curves ####
  
  # Set standard plotting parameters
  # plot_cols <- c( BIC1 = "#37abc8ff", BIC2 = "#ff2a7fff", random_classifier = "grey50",BIC1_boot1 = 'green',BIC1_boot2 = 'red',BIC_3 ='black') # Set plotting colours
  # plot_labels <- c(BIC1 = "BIC1", BIC2 = "BIC2", random_classifier = "Random",
  #                  BIC1_boot1 = 'boot_1' ,BIC1_boot2 = 'boot_2',BIC_3 = 'boot_mean')
  
  # Define overall plotting theme
  theme_set( theme( line = element_line( linewidth = 0.25 ),
                    panel.border = element_rect(colour = "black", linewidth = 0.25, fill = NA), panel.background = element_blank(), panel.grid = element_blank(),
                    axis.text = element_text(size = 6), axis.title = element_text(size = 7),
                    strip.background = element_blank(), strip.text = element_text(size = 7), 
                    legend.title = element_blank(), 
                    legend.background = element_blank(), legend.box.background = element_blank(), legend.text = element_text(size = 6), legend.key.width = unit(4, "mm"), 
                    legend.key = element_rect(fill = NA, colour = NA),
                    plot.title = element_text( size = 7, hjust = 0.5, vjust = 2)))
  
  # Define a function to create a PR curve plot with the AUPRC shown as barchart as an inset
  PR_plot <- function(PR_stats, PR_data, plot_title){
    
    # Plot the AUPRC bar-chart inset as a grob for custom_annotation
    inset <- ggplotGrob(ggplot(PR_stats, aes(x = metric, y = AUPRC, fill = metric))+
                          geom_bar(stat = "identity")+
                          scale_y_continuous( breaks = seq(0,1,0.2))+
                          # scale_fill_manual( values = plot_cols )+
                          theme( legend.position = "none", axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(),
                                 axis.title.y = element_text(size = 6, angle = 90, vjust = 1), axis.text.y = element_text(size = 5)))
    
    # Create the actual plot
    p <- ggplot(PR_data, aes(x = Recall, y = Precision, colour = metric))+
      geom_line( linewidth = 0.25 )+
      ggtitle(plot_title)+
      # scale_colour_manual( values = plot_cols, labels = plot_labels )+
      scale_x_continuous( limits = c(0,1), breaks = seq(0,1,0.2))+
      scale_y_continuous( limits = c(0,1), breaks = seq(0,1,0.2))+
      # theme(legend.position = "none")+
      annotation_custom(grob = inset, xmin = 0.55, xmax = 1, ymin = 0.55, ymax = 1)
    
    return(p)
  }
  
  # Create the plot using my custom plotting function
  PR_stats[,partition:= stringr::str_extract(metric,'.$')]
  PR_data[,partition:= stringr::str_extract(metric,'.$')]
  name_plot = glue::glue('{standard} partition')
  p <- PR_plot(PR_stats = PR_stats, 
               PR_data = PR_data, plot_title = name_plot)
  
  plots_list[[name_plot]]= p
  
}

theme_set(
  theme_bw() +
    theme(legend.position = "top")
)

hyperparams = complete_AUCs[order(-AUPRC)][,head(.SD,1), by = partition][,.(partition,metric)]
hyper_paramlist = list()
AUCs = data.table()
for(i in 1:length(partition_list)){
  print(i)
  DT_partition =fread(here::here('out','datasets', glue::glue('hyperparam_partition__bicor_{i}_surrogate_serule_minplit_cps_balanced.gz')))
  #### Annotate with gold standard ####
  
  # Simplify protein IDs for merging with gold standard
  DT_partition <- DT_partition[ Protein_1 != Protein_2 ]                  # Removing self-associations that arose between isoforms (can't test those anyway)
  DT_partition[ Protein_1 <= Protein_2 ]                        # Validate that A > B still fullfilled   
  plots_list = list()
  for(standard in c('train','test')){
  
    print(standard)
    
    GS_tmp = GS[Class == standard][,set:=NULL]
    setnames(GS_tmp,c('Protein_1','Protein_2','Class'))
    
    GS_tmp[ Protein_1 <= Protein_2, ]     # All pairs sorted A > B, so ready to merge with DT 
    
    # Annotate correlations with gold standards by merging
    DT <- merge(DT_partition, GS, by = c("Protein_1", "Protein_2"))
    DT[, .N, Class]  # Currently about 3% TP pairs
    
    # Downsample FPs to 90% - this is essential to make curves comparable across project
    # n_FP_to_remove <- DT[ Class == "FP", .N ] - DT[ Class == "TP", .N ] * 9
    # DT[ sample( which( Class == "FP" ), n_FP_to_remove ), Class := NA ]    # From the rows where the class is FP, randomly sample a certain subset and turn to NA
    DT <- DT[ !is.na(Class)  ][Class !='']                                             # Remove those rows
    DT[, .N, Class]  # Now 10% TP pairs
    
    # Add a randomised classifier to the input table
    sample_from = DT[,3] |> unlist()
    DT[, random_classifier := sample(sample_from)]
    
    
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
      tmp_curve_dt <- rbind( tmp_curve_dt[ Recall  < 0.05                ][ sample(min(.N, 180))  ],   # Downsample the curves for faster plotting. Keep more points for the low recall areas because they are less robust.
                             tmp_curve_dt[ Recall >= 0.05 & Recall < 0.1 ][ sample(min(.N, 200))  ],
                             tmp_curve_dt[ Recall >= 0.10 & Recall < 0.2 ][sample(min(.N, 200))  ],
                             tmp_curve_dt[ Recall >= 0.20 & Recall < 0.4 ][ sample(min(.N, 200))  ],
                             tmp_curve_dt[ Recall >= 0.40                ][sample(min(.N, 200))  ])
      output <- list( 
        "stats" = data.table( "AUPRC" = tmp_curve$auc.integral, "AUPRC_random" = tmp_curve$rand$auc.integral, "N_TP" = length(tmp_TP), gold_standard = gs, metric = metric), 
        "curve" = tmp_curve_dt)
      return(output)
    }
    
    # Perform the PR analyses using my custom function. 
    PR_curves = list()
    PR_stats = data.table()
    PR_data =data.table()
    for(name in colnames(DT)[-c(1:2,ncol(DT)-1,ncol(DT)-2)]){
      print(name)
      PR_curves[[name]] = PR_curve( DT, "Class", name)
      PR_stats = rbind(PR_stats,
                       PR_curves[[name]][['stats']])
      PR_data = rbind(PR_data,
                      PR_curves[[name]][['curve']])
    }
    PR_stats[,partition:=i]
    PR_stats[,set:=standard]
    AUCs = rbind(AUCs, PR_stats)
    #### Plot the PR curves ####
    
    # Set standard plotting parameters
    # plot_cols <- c( BIC1 = "#37abc8ff", BIC2 = "#ff2a7fff", random_classifier = "grey50",BIC1_boot1 = 'green',BIC1_boot2 = 'red',BIC_3 ='black') # Set plotting colours
    # plot_labels <- c(BIC1 = "BIC1", BIC2 = "BIC2", random_classifier = "Random",
    #                  BIC1_boot1 = 'boot_1' ,BIC1_boot2 = 'boot_2',BIC_3 = 'boot_mean')
    
    # Define overall plotting theme
    theme_set( theme( line = element_line( linewidth = 0.25 ),
                      panel.border = element_rect(colour = "black", linewidth = 0.25, fill = NA), panel.background = element_blank(), panel.grid = element_blank(),
                      axis.text = element_text(size = 6), axis.title = element_text(size = 7),
                      strip.background = element_blank(), strip.text = element_text(size = 7), 
                      legend.title = element_blank(), 
                      legend.background = element_blank(), legend.box.background = element_blank(), legend.text = element_text(size = 6), legend.key.width = unit(4, "mm"), 
                      legend.key = element_rect(fill = NA, colour = NA),
                      plot.title = element_text( size = 7, hjust = 0.5, vjust = 2)))
    
    # Define a function to create a PR curve plot with the AUPRC shown as barchart as an inset
    PR_plot <- function(PR_stats, PR_data, plot_title){
      
      # Plot the AUPRC bar-chart inset as a grob for custom_annotation
      inset <- ggplotGrob(ggplot(PR_stats, aes(x = metric, y = AUPRC, fill = metric))+
                            geom_bar(stat = "identity")+
                            scale_y_continuous( breaks = seq(0,1,0.2))+
                            # scale_fill_manual( values = plot_cols )+
                            theme( legend.position = "none", axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(),
                                   axis.title.y = element_text(size = 6, angle = 90, vjust = 1), axis.text.y = element_text(size = 5)))
      
      # Create the actual plot
      p <- ggplot(PR_data, aes(x = Recall, y = Precision, colour = change))+
        geom_line( linewidth = 0.25 )+
        ggtitle(plot_title)+
        # scale_colour_manual( values = plot_cols, labels = plot_labels )+
        scale_x_continuous( limits = c(0,1), breaks = seq(0,1,0.2))+
        scale_y_continuous( limits = c(0,1), breaks = seq(0,1,0.2))+
        # theme(legend.position = "none")+
        annotation_custom(grob = inset, xmin = 0.55, xmax = 1, ymin = 0.55, ymax = 1)
      
      return(p)
    }
    
    # Create the plot using my custom plotting function
    PR_stats[,partition:= stringr::str_extract(metric,'.$')]
    PR_data[,partition:= stringr::str_extract(metric,'.$')]
    PR_data[,metric_change:=stringr::str_remove(metric,'.$')]
    PR_data[,metric_change:= dplyr::case_when(
      str_detect(metric,'random',negate =F)~ 'random',
      str_detect(metric,'surrogate_0',negate =T)~ 'surrogate',
      str_detect(metric,'splits_20',negate =T)~ 'minsplits',
      str_detect(metric,'cp_0.015',negate =T)~ 'cp',
      
      TRUE~'default'
    ), by = metric]
    PR_data[,change:= dplyr::case_when(
      str_detect(metric,'random',negate =F)~ 'random',
      metric_change == 'surrogate'~str_extract(metric,'surrogate_[:print:]*?_'),
      metric_change == 'minsplits'~str_extract(metric,'splits_[:print:]*?_'),
      metric_change == 'cp'~str_extract(metric,'cp_[:print:]*?_'),
      
      TRUE~'default'
    )]
    
    name_plot = glue::glue('{standard} partition')
    p <- PR_plot(PR_stats = PR_stats,
                 PR_data = PR_data, plot_title = name_plot)+
      facet_grid(partition~metric_change)
    plots_list[[name_plot]]= p
    # name_plot = glue::glue('{gold_standard} partition facet')
    # p <- PR_plot(PR_stats = PR_stats, 
    #              PR_data = PR_data, plot_title = name_plot)+
    #   facet_wrap('partition')
    
    plots_list[[name_plot]]= p
    
  }
}
AUCs[,.(AUPRC,metric, partition,   set)] |> 
  dcast(partition+metric~set, value.var = 'AUPRC') |> 
  ggplot(aes(x = train,y = test, colour = as.character(partition)))+
  geom_point(fill =0.5)+ facet_wrap('partition')+
  geom_abline(slope = 1, intercept = 0)+
  ggtitle('no overfitting in treeclust partition gridsearch')

ggsave(here::here('out','plots','treeclust_partition_train_test.pdf'))
AUCs[,`:=`(splits = str_extract(metric,'splits_[:print:]*?_'),
           cp = str_extract(metric,'cp_[:print:]*?_'),
           serule = str_extract(metric,'serule_[:print:]*?_'),
           surrogate = str_extract(metric,'surrogate_[:print:]*?_'))]
AUCs |> ggplot(aes(x =cp, colour = as.character(serule), y = AUPRC))+
  geom_point()+
  facet_wrap(splits~surrogate)

# making PRcurbes based on the final top3_Scaled score and the original Gold standards
DT_partition = Reduce(function(...) merge(...,by = c("Protein_1", "Protein_2")),list(Whole_ProHD2,
                                                                                     scaled_df_top3))
DT_partition <- DT_partition[ Protein_1 != Protein_2 ]                  # Removing self-associations that arose between isoforms (can't test those anyway)
DT_partition[ Protein_1 <= Protein_2 ]                        # Validate that A > B still fullfilled   
plots_list = list()
GS_tmp = GS[][,set:=NULL]
  setnames(GS_tmp,c('Protein_1','Protein_2','Class'))
  
  GS_tmp[ Protein_1 <= Protein_2, ]     # All pairs sorted A > B, so ready to merge with DT 
  
  # Annotate correlations with gold standards by merging
  DT <- merge(DT_partition, GS_tmp, by = c("Protein_1", "Protein_2"))
  DT[, .N, Class]  # Currently about 3% TP pairs
  
  # Downsample FPs to 90% - this is essential to make curves comparable across project
  # n_FP_to_remove <- DT[ Class == "FP", .N ] - DT[ Class == "TP", .N ] * 9
  # DT[ sample( which( Class == "FP" ), n_FP_to_remove ), Class := NA ]    # From the rows where the class is FP, randomly sample a certain subset and turn to NA
  DT <- DT[ !is.na(Class)  ][Class !='']                                             # Remove those rows
  DT[, .N, Class]  # Now 10% TP pairs
  
  # Add a randomised classifier to the input table
  sample_from = DT[,3] |> unlist()
  DT[, random_classifier := sample(sample_from)]
  
  
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
    tmp_curve_dt <- rbind( tmp_curve_dt[ Recall  < 0.05                ][ sample(min(.N, 180)) ],   # Downsample the curves for faster plotting. Keep more points for the low recall areas because they are less robust.
                           tmp_curve_dt[ Recall >= 0.05 & Recall < 0.1 ][ sample(min(.N, 200)) ],
                           tmp_curve_dt[ Recall >= 0.10 & Recall < 0.2 ][ sample(min(.N, 200))  ],
                           tmp_curve_dt[ Recall >= 0.20 & Recall < 0.4 ][ sample(min(.N, 200))  ],
                           tmp_curve_dt[ Recall >= 0.40                ][sample(min(.N, 200)) ])
    output <- list( 
      "stats" = data.table( "AUPRC" = tmp_curve$auc.integral, "AUPRC_random" = tmp_curve$rand$auc.integral, "N_TP" = length(tmp_TP), gold_standard = gs, metric = metric), 
      "curve" = tmp_curve_dt)
    return(output)
  }
  
  # Perform the PR analyses using my custom function. 
  PR_curves = list()
  PR_stats = data.table()
  PR_data =data.table()
  for(name in colnames(DT)[-c(1:2,ncol(DT)-1)]){
    print(name)
    PR_curves[[name]] = PR_curve( DT, "Class", name)
    PR_stats = rbind(PR_stats,
                     PR_curves[[name]][['stats']])
    PR_data = rbind(PR_data,
                    PR_curves[[name]][['curve']])
  }
  theme_set( theme( line = element_line( linewidth = 0.25 ),
                    panel.border = element_rect(colour = "black", linewidth = 0.25, fill = NA), panel.background = element_blank(), panel.grid = element_blank(),
                    axis.text = element_text(size = 6), axis.title = element_text(size = 7),
                    strip.background = element_blank(), strip.text = element_text(size = 7), 
                    legend.title = element_blank(), 
                    legend.background = element_blank(), legend.box.background = element_blank(), legend.text = element_text(size = 6), legend.key.width = unit(4, "mm"), 
                    legend.key = element_rect(fill = NA, colour = NA),
                    plot.title = element_text( size = 7, hjust = 0.5, vjust = 2)))
  
  # Define a function to create a PR curve plot with the AUPRC shown as barchart as an inset
  PR_plot <- function(PR_stats, PR_data, plot_title){
    
    # Plot the AUPRC bar-chart inset as a grob for custom_annotation
    inset <- ggplotGrob(ggplot(PR_stats, aes(x = metric, y = AUPRC, fill = metric))+
                          geom_bar(stat = "identity")+
                          scale_y_continuous( breaks = seq(0,1,0.2))+
                          # scale_fill_manual( values = plot_cols )+
                          theme( legend.position = "none", axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(),
                                 axis.title.y = element_text(size = 6, angle = 90, vjust = 1), axis.text.y = element_text(size = 5)))
    
    # Create the actual plot
    p <- ggplot(PR_data, aes(x = Recall, y = Precision, colour = metric))+
      geom_line( linewidth = 0.25 )+
      ggtitle(plot_title)+
      # scale_colour_manual( values = plot_cols, labels = plot_labels )+
      scale_x_continuous( limits = c(0,1), breaks = seq(0,1,0.2))+
      scale_y_continuous( limits = c(0,1), breaks = seq(0,1,0.2))+
      # theme(legend.position = "none")+
      annotation_custom(grob = inset, xmin = 0.55, xmax = 1, ymin = 0.55, ymax = 1)
    
    return(p)
  }
  
  
  p <- PR_plot(PR_stats = PR_stats,
               PR_data = PR_data, plot_title = 'various covariations')+
    theme(legend.position = 'bottom')
  
ggsave(here::here('out','plots','whole_BIC_PCC_TC.pdf'),p)
