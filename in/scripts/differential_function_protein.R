# changes to function of single protein
# BiocManager::install("GEOquery")
library(diffcoexp); library(ggpmisc); library(data.table)
library(UpSetR); library(org.Calbicans.eg.db)
library(clusterProfiler); library(future);library(stringr)
DT_normalised = fread(here::here('out',
                                 'datasets',
                                 'covariation_stressed_treatments_normalised.gz'))

DT_normalised[bicor>0.57, .N ,by = Condition]
DT_normalised[,Pair := paste(Protein_1,Protein_2,sep ='_')]


significant_pairs = DT_normalised[bicor>0.785]
significant_pairs[,.N,by = Condition]
list_significant_pairs = list()
for(i in unique(DT_normalised$Condition)){
  print(i)
  list_significant_pairs[[i]] <- DT_normalised[bicor>0.57 & Condition==i & Pair %in% significant_pairs$Pair ,Pair]

}
UpSetR::upset(UpSetR::fromList(list_significant_pairs), nsets = 11,
              order.by = 'freq')
grid::grid.text("high confidence pairs (bicor >0.785) in at least 1, \npresent in >0.57 in all",x = 0.65, y=0.95)

significant_pairs = DT_normalised[between(bicor,0.57,0.785)]
significant_pairs[,.N,by = Condition]
list_significant_pairs = list()
for(i in unique(DT_normalised$Condition)){
  print(i)
  list_significant_pairs[[i]] <- DT_normalised[bicor>0.57 & Condition==i & Pair %in% significant_pairs$Pair ,Pair]
  
}
UpSetR::upset(UpSetR::fromList(list_significant_pairs), nsets = 11,
              order.by = 'freq')
grid::grid.text("high confidence pairs (0.57<bicor <0.785) in at least 1, \npresent in >0.57 uniquely",x = 0.65, y=0.95)


DT_pair = DT_normalised[Condition %in% c('SC_flc','SM_lys') & bicor>0.57,.(Protein_1,Protein_2)] |> 
  unique()
DT_pair= merge(DT_pair,
               DT_normalised[Condition %in% c('SC_flc','SM_lys')],
               by = c('Protein_1','Protein_2'))
DT_pair = DT_pair |> dcast(Protein_1+Protein_2~Condition,value.var = 'bicor')
DT_pair[,differential_corr := SC_flc-SM_lys]
DT_pair_rev = copy(DT_pair)
DT_pair_rev[,`:=`(Protein_1 = Protein_2,
                  Protein_2 = Protein_1)]
DT_pair = rbind(DT_pair,
                DT_pair_rev)
diff_functions = DT_pair[,.(mean_abs = mean(abs(differential_corr),na.rm = T),
                            N_pairs = .N), by = Protein_1] 
diff_functions |> 
  ggplot(aes(y = N_pairs, x = mean_abs, label = Protein_1))+
  geom_point()+
  ggrepel::geom_text_repel(data = diff_functions[N_pairs>5 & mean_abs>0.5])

examples = c('A0A1D8PEH5','A0A1D8PRM4')
examples_interactors = DT_normalised[(Protein_1 == examples[1] |
                                       Protein_2 == examples[1]) &
                                       bicor>0.4 &
                                       Condition %in% c('SC_flc','SM_lys')]
examples_interactors[,Partner:= fifelse(Protein_1 == examples[1] ,
                                        Protein_2, Protein_1)]
UpSetR::upset(UpSetR::fromList(list(SC_flc=examples_interactors[Condition == 'SC_flc',Partner],
                                    SM_lys = examples_interactors[Condition == 'SM_lys',Partner])) )

examples_interactors = DT_normalised[(Protein_1 == examples[2] |
                                        Protein_2 == examples[2]) &
                                       bicor>0.4 &
                                       Condition %in% c('SC_flc','SM_lys')]
examples_interactors[,Partner:= fifelse(Protein_1 == examples[2] ,
                                        Protein_2, Protein_1)]
UpSetR::upset(UpSetR::fromList(list(SC_flc=examples_interactors[Condition == 'SC_flc',Partner],
                                    SM_lys = examples_interactors[Condition == 'SM_lys',Partner])) )


# plotting_protein_specific_correlation

DT_plot= DT_normalised[Condition %in% c('SC_flc','SM_lys')]
DT_plot = DT_plot |> dcast(Protein_1+Protein_2~Condition,value.var = 'bicor')
DT_pair_rev = copy(DT_plot)
DT_pair_rev[,`:=`(Protein_1 = Protein_2,
                  Protein_2 = Protein_1)]
DT_plot = rbind(DT_plot,
                DT_pair_rev)
DT_plot[Protein_1 == 'A0A1D8PJ34'] |> ggplot(aes(x = SC_flc, y = SM_lys ))+
  geom_point()


# diffcoexp

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
conditions = colnames(df)  |> str_extract('_S[:print:]*') |> str_remove('^_') |> unique()

list_dataset = list()
for(i in conditions){
  print(i)
  dt_tmp = df[,str_detect(colnames(df),paste0(i,'$'))]
  # Transpose and median normalise the rows to avoid spurious correlations
  tmp_medians <- apply( dt_tmp, 1, median, na.rm = TRUE )  
  dt_tmp <- sweep( dt_tmp, 1, tmp_medians, FUN = "-" )
  dt_tmp = dt_tmp[matrixStats::rowSums2(is.na(dt_tmp))<20,]
  list_dataset[[i]]= dt_tmp
  
}
list_dataset |> purrr::map_dbl(ncol)
list_dataset |> purrr::map_dbl(nrow)

library(diffcoexp)
allowWGCNAThreads(nThreads = 40)
common_genes = intersect(rownames(list_dataset$SC_flc),
          rownames(list_dataset$SM_lys))


res=diffcoexp(exprs.1 = list_dataset$SC_flc[common_genes,],
              q.method  ='bonferroni',
              rth = 0.7,
              qth = 0.0001,  q.diffth = 0.0001, 
              q.dcgth = 0.0001,
              r.diffth = 0.3,
              # r.diffth = 
              exprs.2 = list_dataset$SM_ly[common_genes,], r.method = "pearson" )
res$DCGs |> View()

# plotting_protein_specific_correlation

# DT_normalised = fread(here::here('out',
#                                  'datasets',
#                                  'covariation_stressed_treatments_normalised.gz'))
DT_plot= DT_normalised[Condition %in% c('SC_flc','SM_lys')]
DT_plot = DT_plot |> dcast(Protein_1+Protein_2~Condition,value.var = 'bicor')
DT_pair_rev = copy(DT_plot)
DT_pair_rev[,`:=`(Protein_1 = Protein_2,
                  Protein_2 = Protein_1)]
DT_plot = rbind(DT_plot,
                DT_pair_rev)
Protein_of_interest = 'Q5AD22'
plot_of_interest = DT_plot[Protein_1 == Protein_of_interest]
partners = res$DCLs |> as.data.table()
partners = partners[Gene.1 ==Protein_of_interest | Gene.2 ==Protein_of_interest]
partners[,Protein_2:= fifelse(Gene.1 ==Protein_of_interest, Gene.2, Gene.1 )]
partners = merge(plot_of_interest,
      partners, by = 'Protein_2', all = T) 
  ggplot(partners,aes(x = SM_lys , y =  SC_flc, colour = type ))+
  geom_point(alpha = 0.1)+
  geom_point(data = partners[!is.na(type)], alpha = 1)+
    lims(x  = c(-1,1),y = c(-1,1))

  # weighted jaccard score
  ## simulate data
  # nr <- 2291; nc <- 265
  # set.seed(420)
  # input_m <- matrix(rnorm(nr * nc), nrow = nr, ncol = nc)
  # input_m[1:5, 1:5]
  # input_m[input_m<0.3] <-  0
  # 
  jaccardMinem <- function(input_m) {
    require(data.table)
    require(matrixStats)
    
    samples <- 1:ncol(input_m)
    comb <- CJ(samples, samples)
    comb[, i := .I]
    comb <- melt(comb, 'i')
    setorder(comb, value)
    v2 <- paste0("V", 1:2)
    comb[, variable2 := v2 , keyby = i]
    comb2 <- dcast(comb, i ~ variable2, value.var = 'value')
    combUnique <- unique(comb2, by = c('V1', 'V2'))
    
    XX <- apply(combUnique[, -'i'], 1, function(x) {
      x2 <- rowRanges(input_m, cols = x)
      s <- colSums2(x2)
      s[1] / s[2]
    })
    
    set(combUnique, j = 'xx', value = XX)
    rez2 <- merge(comb2, combUnique[, -'i'], by = c('V1', 'V2'), all.x = T)
    setorder(rez2, i)
    rez2 <- array(rez2$xx, dim = rep(ncol(input_m), 2))
    rownames(rez2) <- colnames(input_m)
    colnames(rez2) <- colnames(input_m)
    rez2
  }
  # jacc_interactor_similarity = jaccardMinem(input_m)
  # jacc_interactor_similarity |> dim()

  # jacc_interactor_similarity |> pheatmap::pheatmap()
  
  # matrix per protein across conditions
  signficant_one_conditions = DT_normalised[bicor>0.57]
  all_pairs_threshold = rbind(DT_normalised[bicor>0.3,.(Protein_1,Protein_2,Condition,bicor)],
                              DT_normalised[,.(bicor = mean(bicor,na.rm= T)), by = .(Protein_1,Protein_2)][,Condition := 'avg_across'][bicor>0.3])
  all_proteins_across_conditions = data.table()
  for(protein in unique(c(signficant_one_conditions$Protein_1,signficant_one_conditions$Protein_2)) ){
    print(protein)
    interactors = signficant_one_conditions[Protein_1 == protein|
                                              Protein_2 == protein]
    interactors_avg_bicor = interactors[,.(avg_bicor = mean(bicor),
                   total_partners = .N)]
    
    interactors= interactors[, .(Protein_1,Protein_2)] |> unique()
   matrix_partners =  merge(all_pairs_threshold,
          interactors, by = c('Protein_1','Protein_2')) |> 
      dcast(Protein_1 + Protein_2 ~Condition, value.var = 'bicor')
   matrix_partners = matrix_partners[,-c(1:2)] |> 
     as.matrix()
   matrix_partners[is.na(matrix_partners)] <- 0
   tmp = jaccardMinem(matrix_partners)
   # pheatmap::pheatmap(tmp)
   tmp <- as.data.table( reshape2::melt( as.matrix( tmp  )))                                   # Turn distance matrix into a pair-wise data.table
   tmp <- tmp[, .(Condition_1 = as.character(Var1), Condition_2 = as.character(Var2), value ) ]   # Rename and change to character
   tmp <- tmp[ Condition_1 > Condition_2 ][,`:=`(Protein =  protein,
                                                 n_pairs = nrow(interactors),
                                                 avg_interactor_bicor= interactors_avg_bicor$avg_bicor,
                                                 total_significant_pairs = interactors_avg_bicor$total_partners)]                                                       # Remove redundant pairs by keeping only A > B, removing A == B and B < A pairs
   all_proteins_across_conditions = rbind(all_proteins_across_conditions,
                                          tmp)
  }
  all_proteins_across_conditions$value |> hist()
  all_proteins_across_conditions_avg = all_proteins_across_conditions[,.(median_wjcc = median(value),
                                    sd_wjcc = sd(value),
                                    avg_interactor_bicor = mean(avg_interactor_bicor ),
                                    total_significant_pairs = mean(total_significant_pairs),
                                    N_partners = median(n_pairs)), by = Protein]
    ggplot(all_proteins_across_conditions_avg[N_partners>10],aes(x = median_wjcc, colour =log2(total_significant_pairs),
                                                                 y =sd_wjcc, label = Protein, alpha = 0.5 ))+
    geom_point(aes(size = N_partners))+
      theme_bw()+
    ggrepel::geom_text_repel(data =all_proteins_across_conditions_avg[!between(sd_wjcc,0.1,0.2) &  
                                                                        total_significant_pairs>20] )+
      scico::scale_color_scico(palette = 'batlow')
  
    # examples 
    # protein = 'A0A1D8PFE8'
    # protein = 'A0A1D8PU83'
    protein = 'A0A1D8PGN9'
    # protein = 'Q59S50'
    # protein = 'Q59WJ3'
    interactors = signficant_one_conditions[Protein_1 == protein|
                                              Protein_2 == protein, .(Protein_1,Protein_2)] |> 
      unique()
    matrix_partners =  merge(all_pairs_threshold,
                             interactors, by = c('Protein_1','Protein_2')) |> 
      dcast(Protein_1 + Protein_2 ~Condition, value.var = 'bicor')
    matrix_partners = matrix_partners[,-c(1:2)] |> 
      as.matrix()
    matrix_partners[is.na(matrix_partners)] <- 0
    tmp = jaccardMinem(matrix_partners)
    pheatmap::pheatmap(tmp, 
                       color = scico::scico(10, palette = 'batlow' ),
                       breaks = seq(0,1,0.1),
                       main = glue::glue('{protein} interactor agreement across conditions'))
    all_interactors = merge(DT_normalised,
          interactors, by = c('Protein_1','Protein_2'))
    all_interactors[,interactor := fifelse(Protein_1 == protein,
                                       Protein_2, Protein_1)]
    all_interactors[,avg_bicor := mean(bicor), by = interactor ]
    all_interactors[Condition != 'avg_bicor' & bicor >0.57,
                    N_conditions_significant := 1 ]
    all_interactors[,N_conditions_significant := sum(N_conditions_significant,na.rm =T), by = interactor]
    all_interactors |> ggplot(aes(x = avg_bicor , y = bicor, colour = N_conditions_significant))+
      geom_point(alpha = 0.9)+
      facet_wrap('Condition')+
      lims(x = c(-1,1), y = c(-1,1))+
      # stat_poly_line() +
      stat_poly_eq() +
      theme_bw()+
      geom_abline(intercept = 0, slope = 1, linetype="dotted")+
      annotate("rect", xmin = 0.57, xmax = 1, ymin = -1, ymax = 1,
               alpha = .2, colour = 'darkred')+
      annotate("rect", xmin = -1, xmax = 1, ymin = 0.57, ymax = 1,
               alpha = .2, colour = 'darkred')+
      scico::scale_color_scico(palette = 'batlow')+
      ggtitle(glue::glue('{protein} interactor agreement across conditions'))
   
    overlap_interactors = all_interactors[Condition != 'avg_bicor' & bicor >0.57,.(Condition,interactor)] |> 
      as.data.frame()

    UpSetR::upset(UpSetR::fromList(    split(overlap_interactors$interactor,overlap_interactors$Condition, drop = T)),
                  nsets = 20)
    
    all_interactors[,.(max_cor = max(bicor)), by = interactor]
    
    all_proteins_across_conditions[Condition_1 =='avg_across'    |
                                     Condition_2 =='avg_across'  ][n_pairs>10] |> View()
    
    conditions =  all_proteins_across_conditions[Condition_1 =='avg_across'    |
                                                   Condition_2 =='avg_across' 
                                                 ][,Condition := fifelse(Condition_1 =='avg_across',
                                                                         Condition_2,
                                                                         Condition_1)]
    # Uniprot mapping
    Uniprot_annot = fread(here::here('in','datasets','uniprotkb_taxonomy_id_237561_2025_07_08.tsv.gz'))
    Uniprot_annot =Uniprot_annot[GeneID != '']
    Uniprot_annot[,CGD:=str_remove(CGD,';[:print:]*$')]
    Uniprot_annot[,GeneID:=str_remove(GeneID,';[:print:]*$')]
    
    condition_specific_conditions = data.table()
    for(i in unique(conditions$Condition)){
      print(i)
      conditions_tmp = conditions[Condition ==i][n_pairs >10,.(Protein,value)]
      conditions_tmp= merge(conditions_tmp,Uniprot_annot[,.(Entry,GeneID)],
            by.x ='Protein', by.y = 'Entry')
      
      conditions_tmp= conditions_tmp[(order(value)),.( GeneID,value)] |> 
        tibble::deframe()
      conditions_tmp = scale(conditions_tmp)* (-1) |> 
        tibble::deframe()
      agreement_avg = conditions_tmp |> as.vector()
      names(agreement_avg) = rownames(conditions_tmp)
      
      ego3 <- gseGO(geneList     = agreement_avg,
                    OrgDb        = org.Calbicans.eg.db,
                    ont          = "ALL",
                    minGSSize    = 20,
                    maxGSSize    = 500,
                    pvalueCutoff = 0.05)
      results = ego3@result |> as.data.table()
      condition_specific_conditions= rbind(condition_specific_conditions,
                                           results[,Condition := i])
      # conditions_tmp[total_significant_pairs  >10] |> ggplot(aes(x = n_pairs, y= value))+
      #     geom_point(alpha = 0.1)
      
      
    }
    condition_specific_conditions_plot = condition_specific_conditions[p.adjust<0.01] 
      ggplot(condition_specific_conditions_plot, aes(y = NES, x = Condition, label = Description))+
      geom_jitter(height = 0.1, width = 0)+
        theme_bw()+
        ggrepel::geom_text_repel(data  = condition_specific_conditions_plot[order(-NES)][,head(.SD,2), by = Condition], 
                                 max.overlaps = 10)+
        ggrepel::geom_text_repel(data  = condition_specific_conditions_plot[order(NES)][,head(.SD,2), by = Condition], 
                               max.overlaps = 5)+
        ggtitle('Condition specific changes in covariation partner agreement')
      
      # comparing oxidatitve phosphorylation across conditions
      go_OI = 'proteasomal protein catabolic process'
      proteins_associations = condition_specific_conditions_plot[Description == go_OI,core_enrichment ] |> 
        str_split('/') |> unlist() |> unique()
      proteins_associations = Uniprot_annot[GeneID %in% proteins_associations,Entry]
      terms_conditions_specific = DT_normalised[][Protein_1 %in% proteins_associations &
                        Protein_2 %in% proteins_associations, type:=go_OI]
      terms_conditions_specific[is.na(type),type:= 'other']
      terms_conditions_specific = terms_conditions_specific[,.SD[sample(.N, min(.N,10000))], by = .(Condition,type)]
      terms_conditions_specific |> ggplot(aes(x = Condition, y = bicor, colour = type))+
        geom_boxplot()+
        theme_bw()+
        facet_wrap('type')
      
      
      
# change function protein
      # according to top interactors
      # Uniprot mapping
      Uniprot_annot = fread(here::here('in','datasets','uniprotkb_taxonomy_id_237561_2025_07_08.tsv.gz'))
      Uniprot_annot[,CGD:=str_remove(CGD,';[:print:]*$')]
      Uniprot_annot[,GeneID:=str_remove(GeneID,';[:print:]*$')]
    universe_geneID = Uniprot_annot[GeneID != '',.(Entry,GeneID)]
    
    # finding the top 100 interactors per protein
    DT_normalised_rev = copy(DT_normalised)
    DT_normalised_rev[,`:=`(Protein_1 = Protein_2,
                            Protein_2 = Protein_1)]
    DT_normalised_rev = rbind(DT_normalised_rev,
                              DT_normalised)
    covariation_partners =DT_normalised_rev[order(-bicor)
                      ][,.(Protein_1,Protein_2,bicor,Condition)
                        ][!is.na(bicor),head(.SD,50), by = .(Protein_1, Condition)]
    setkey(covariation_partners,'Protein_1')
    
    
    calc_pairwise_overlaps <- function(sets) {
      # Ensure that all sets are unique character vectors
      sets_are_vectors <- vapply(sets, is.vector, logical(1))
      if (any(!sets_are_vectors)) {
        stop("Sets must be vectors")
      }
      sets_are_atomic <- vapply(sets, is.atomic, logical(1))
      if (any(!sets_are_atomic)) {
        stop("Sets must be atomic vectors, i.e. not lists")
      }
      sets <- lapply(sets, as.character)
      is_unique <- function(x) length(unique(x)) == length(x)
      sets_are_unique <- vapply(sets, is_unique, logical(1))
      if (any(!sets_are_unique)) {
        stop("Sets must be unique, i.e. no duplicated elements")
      }
      
      n_sets <- length(sets)
      set_names <- names(sets)
      n_overlaps <- choose(n = n_sets, k = 2)
      
      vec_name1 <- character(length = n_overlaps)
      vec_name2 <- character(length = n_overlaps)
      vec_num_shared <- integer(length = n_overlaps)
      vec_overlap <- numeric(length = n_overlaps)
      vec_jaccard <- numeric(length = n_overlaps)
      overlaps_index <- 1
      
      for (i in seq_len(n_sets - 1)) {
        name1 <- set_names[i]
        set1 <- sets[[i]]
        for (j in seq(i + 1, n_sets)) {
          name2 <- set_names[j]
          set2 <- sets[[j]]
          
          set_intersect <- set1[match(set2, set1, 0L)]
          set_union <- .Internal(unique(c(set1, set2), incomparables = FALSE,
                                        fromLast = FALSE, nmax = NA))
          num_shared <- length(set_intersect)
          overlap <- num_shared / min(length(set1), length(set2))
          jaccard <- num_shared / length(set_union)
          
          vec_name1[overlaps_index] <- name1
          vec_name2[overlaps_index] <- name2
          vec_num_shared[overlaps_index] <- num_shared
          vec_overlap[overlaps_index] <- overlap
          vec_jaccard[overlaps_index] <- jaccard
          
          overlaps_index <- overlaps_index + 1
        }
      }
      
      result <- data.frame(name1 = vec_name1,
                           name2 = vec_name2,
                           # n_members_1 = length_vec_1,
                           # n_members_2 = length_vec_2,
                           num_shared = vec_num_shared,
                           overlap = vec_overlap,
                           jaccard = vec_jaccard,
                           stringsAsFactors = FALSE)
      return(result)
    }
    top_partner_similarity= data.table()
    for(i in unique(covariation_partners$Protein_1)){
      print(i)
     partners_tmp =  covariation_partners[Protein_1 ==i]
     consensus = partners_tmp[,.N, by = Protein_2
                              ][order(-N)][,head(.SD,50)
                                           ][,Protein_2]
     partners = split(partners_tmp$Protein_2,partners_tmp$Condition)
     partners[['consensus']] <- consensus
     
     partners_overlap = calc_pairwise_overlaps(partners) |> 
       as.data.table()
     top_partner_similarity = rbind(top_partner_similarity,
                                    partners_overlap[,Protein:=i])
    }
    statscond = covariation_partners[,.(Min_cor = min(bicor),
                            avg_cor = mean(bicor)),by = .(Protein_1,Condition) ]
    top_partner_similarity[,Condition:= fifelse(name1 =='consensus',name2,name1)] 
    top_partner_similarity_stats =  merge(top_partner_similarity[name1 =='consensus' |
                                                                   name2 == 'consensus'],
                                          statscond,
          by.x = c('Protein','Condition'),by.y = c('Protein_1','Condition'))
         
    
    top_partner_similarity_stats |> 
      ggplot(aes(x = overlap, y = avg_cor))+
     geom_hex()+
      theme_bw()+
      facet_wrap('Condition')+
      labs(x = 'interactor overlap per protein', y = 'mean interactor covariance')+
      ggtitle('The agreement of top interactors between conditions - same function',
              subtitle = 'agrees with the strength of covariation')
      ggsave(here::here('out','plots','overlap_interator_top_conditions.pdf'))
    
      Protein_OI = 'A0A1D8PKD5'
      Protein_OI_dt = DT_normalised[Protein_1 == Protein_OI|
                             Protein_2 == Protein_OI][order(-bicor)] 
    Protein_OI_dt[,Partner:= fifelse(Protein_1 == Protein_OI,Protein_2,Protein_1)]
    SC_edta_interactors = universe_geneID[Entry %in%  
                                            Protein_OI_dt[Condition=='SC_edta'
                                                          ][,head(.SD,50)][,Partner],
                                          GeneID]
    universe_edta = universe_geneID[Entry %in%  
                                      Protein_OI_dt[Condition=='SC_edta'
                                      ][,Partner],
                                    GeneID]
    SM_lys_interactors = universe_geneID[Entry %in%  
                                            Protein_OI_dt[Condition=='SM_lys'
                                            ][,head(.SD,50)][,Partner],
                                          GeneID]
    
    universe_SM_lys = universe_geneID[Entry %in%  
                                        Protein_OI_dt[Condition=='SM_lys'
                                        ][,Partner],
                                      GeneID]
        ggplot(Protein_OI_dt,aes(x = avg_cor, y = bicor))+
        geom_hex()+
        facet_wrap('Condition')+
        theme_bw()+
          geom_abline(slope = 1, intercept = 0)+
        labs(x = 'Avg cov interactor across conditions',
             y = 'partner covariation in specific condition')
      ggsave(here::here('out','plots',glue::glue('{Protein_OI}_overlap_interator_top_conditions.pdf')))
      
      ego_edta <- enrichGO(gene          =SC_edta_interactors,
                      universe      = universe_edta,  
                      OrgDb         = org.Calbicans.eg.db,
                      ont           = "ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      # keyType = 'UNIPROT',
                      readable      = TRUE)
      ego_SM_lys <- enrichGO(gene          =SM_lys_interactors,
                           universe      = universe_SM_lys,  
                           OrgDb         = org.Calbicans.eg.db,
                           ont           = "ALL",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.01,
                           qvalueCutoff  = 0.05,
                           # keyType = 'UNIPROT',
                           readable      = TRUE)
      dotplot(ego_edta)+
        ggtitle(glue::glue('{Protein_OI}_SC_edta_top50 enrichment'))
      ggsave(here::here('out','plots',glue::glue('{Protein_OI}_SC_edta_top50 enrichment.pdf')))
      
      dotplot(ego_SM_lys)+
        ggtitle(glue::glue('{Protein_OI}_SM_lys_top50 enrichment'))
      ggsave(here::here('out','plots',glue::glue('{Protein_OI}_SM_lys_top50 enrichment.pdf')))
      
      covariation_partners = merge(DT_normalised,
                                universe_geneID[,.(Entry,GeneID)],
                                   by.x = 'Protein_1',
                                   by.y = 'Entry') |> 
        merge(universe_geneID[,.(Entry,GeneID)],
              by.x = 'Protein_2',
              by.y = 'Entry') 
      covariation_partners = covariation_partners[GeneID.x != ''
      ][GeneID.y != ''   ]
      enrichment_function <- function(Uniprot,Condition,Universe,covariations){
        i  = universe_geneID[Entry == Uniprot,GeneID]
        interactors = covariations[GeneID.x == i|
                                             GeneID.y == i 
        ][,Interactor:= fifelse(GeneID.x==i,GeneID.y,GeneID.x) ][,Interactor]
        if(length(interactors)>2){
          ego <- enrichGO(gene          =interactors ,
                          universe      = universe_geneID$GeneID,  
                          OrgDb         = org.Calbicans.eg.db,
                          ont           = "ALL",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.01,
                          qvalueCutoff  = 0.05,
                          # keyType = 'UNIPROT',
                          readable      = TRUE)
          if(!(is.null(ego)) ){
            term_OI = ego@result |> 
              subset(p.adjust<0.05) |> 
              dplyr::arrange(-Count) |> as.data.table()
            if(nrow(term_OI)>1){
              term_OI[,Protein:=Uniprot]
              term_OI[,Condition:=Condition]
              fwrite(term_OI, here::here('out','datasets',glue::glue('enrichment_terms_bicor_conditionss.csv')), append = T)
              
            }
            
          }
        }
      }
      
      plan(sequential)
      plan(multisession, workers = 40)
      
      for(i in unique(covariation_partners$Condition)){
        print(i)
      furrr::future_walk(.x = universe_geneID$Entry, ~enrichment_function(Uniprot = .x,
                                                                          Condition = i,
                                                                          Universe = universe_geneID,
                                                                          covariations =covariation_partners[Condition ==i] ))
      }  
      
      
      # change function protein
      # according to signficant interactors
      # Uniprot mapping
      Uniprot_annot = fread(here::here('in','datasets','uniprotkb_taxonomy_id_237561_2025_07_08.tsv.gz'))
      Uniprot_annot[,CGD:=str_remove(CGD,';[:print:]*$')]
      Uniprot_annot[,GeneID:=str_remove(GeneID,';[:print:]*$')]
      covar_quantiles = DT_normalised[!is.na(bicor),bicor] |> quantile(probs = seq(0, 1, 0.001), type = 5)
      
      covariation_partners <- DT_normalised[bicor>=covar_quantiles["99.0%"],.(Protein_1,Protein_2,bicor,Condition)]
      universe = unique(c(covariation_partners$Protein_1,covariation_partners$Protein_2))
      # number of proteins with at least 1 covariation Partner
      universe |> length()
      # universe are all proteins with at least 1 covariation partner
      universe_geneID = Uniprot_annot[Entry %in%universe ,.(Entry,GeneID)]
      Uniprot_annot =Uniprot_annot[GeneID != '']
      
      covariation_partners = merge(covariation_partners,
                                   Uniprot_annot[,.(Entry,GeneID)],
                                   by.x = 'Protein_1',
                                   by.y = 'Entry') |> 
        merge(Uniprot_annot[,.(Entry,GeneID)],
              by.x = 'Protein_2',
              by.y = 'Entry') 
      covariation_partners = covariation_partners[GeneID.x != ''
      ][GeneID.y != ''   ]
      enrichment_function <- function(Uniprot,Condition,Universe,covariations){
        i  = universe_geneID[Entry == Uniprot,GeneID]
        interactors = covariations[GeneID.x == i|
                                     GeneID.y == i 
        ][,Interactor:= fifelse(GeneID.x==i,GeneID.y,GeneID.x) ][,Interactor]
        if(length(interactors)>2){
          ego <- enrichGO(gene          =interactors ,
                          universe      = universe_geneID$GeneID,  
                          OrgDb         = org.Calbicans.eg.db,
                          ont           = "ALL",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.01,
                          qvalueCutoff  = 0.05,
                          # keyType = 'UNIPROT',
                          readable      = TRUE)
          if(!(is.null(ego)) ){
            term_OI = ego@result |> 
              subset(p.adjust<0.05) |> 
              dplyr::arrange(-Count) |> as.data.table()
            if(nrow(term_OI)>1){
              term_OI[,Protein:=Uniprot]
              term_OI[,Condition:=Condition]
              fwrite(term_OI, here::here('out','datasets',glue::glue('enrichment_terms_bicor_conditionss.csv')), append = T)
              
            }
            
          }
        }
      }
      
      plan(sequential)
      plan(multisession, workers = 40)
      
      for(i in unique(covariation_partners$Condition)){
        print(i)
        furrr::future_walk(.x = universe_geneID$Entry, ~enrichment_function(Uniprot = .x,
                                                                            Condition = i,
                                                                            Universe = universe_geneID,
                                                                            covariations =covariation_partners[Condition ==i] ))
      }  
      