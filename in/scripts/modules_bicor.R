# predicting modules
library(clusterProfiler);library(org.Calbicans.eg.db); library(data.table)
bicor_covar = fread(here::here('out','datasets','bicor_all_samples.gz'))
# cluster enrichment
# Uniprot mapping
Uniprot_annot = fread(here::here('in','datasets','uniprotkb_taxonomy_id_237561_2025_07_08.tsv.gz'))
Uniprot_annot[,CGD:=str_remove(CGD,';[:print:]*$')]
Uniprot_annot[,GeneID:=str_remove(GeneID,';[:print:]*$')]
Uniprot_annot[,`Gene Names`:=str_remove(`Gene Names`,' [:print:]*$')]
Uniprot_annot[`Gene Names`=='',`Gene Names` := Entry]
covar_quantiles = bicor_covar[!is.na(corr_bicor),corr_bicor] |> quantile(probs = seq(0, 1, 0.001), type = 5)
covariation_partners <- bicor_covar[corr_bicor>=0.56,.(Protein_1,Protein_2,corr_bicor)]


fwrite(bicor_covar[corr_bicor>0.56], col.names = F,sep = '\t',
       here::here('out','datasets','cluster_one_input.txt'))
command = 'java -jar /home/v1skourt/Downloads/cluster_one-1.0.jar /home/v1skourt/candida/out/datasets/cluster_one_input.txt -output-format csv --min-density 0.4'
clusterONE = system(command, intern = T)

clusters = data.table(cols =clusterONE[-1])
clusterONE = tidyr::separate_wider_delim(clusters,cols = 'cols', delim =',',
                            names = unlist(str_split(clusterONE[1],',',simplify = T))) 
clusterONE =  as.data.table(clusterONE)
char_cols = colnames(clusterONE)[c(2:7)]
clusterONE[,(char_cols) := lapply(.SD, as.numeric), .SDcols = char_cols]

remove_backlash <- function(x){
  gsub('[(\")(\\)]', '', deparse(x))
}
clusterONE[,Members:=remove_backlash(Members),by = Cluster]
fwrite(clusterONE,here::here('out','datasets','all_clusterONE_modules.gz'))
sets_of_complexes = purrr::map(.x = clusterONE$Members,
           ~as.vector(str_split(.x,' ',simplify = T)))
names(sets_of_complexes) <- clusterONE$Cluster
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

clusters_overlap = calc_pairwise_overlaps(sets_of_complexes) |> 
  as.data.table()
clusters_overlap = merge(clusters_overlap,clusterONE[,.(Cluster,Size)], by.x= 'name1', by.y ='Cluster') |> 
  merge(clusterONE[,.(Cluster,Size)], by.x= 'name2', by.y ='Cluster') 
clusters__to_remove = clusters_overlap[jaccard>0.2
                 ][,to_remove:= fifelse(Size.x>Size.y,name2,name1)
                   ][,to_remove] |> unique()
clusterONE[,to_remove := Cluster %in% clusters__to_remove]
clusterONE$to_remove |> table()
clusterONE |> ggplot(aes(y = Size, x = to_remove))+
  geom_boxplot(alpha = 0.5)
universe = unique(c(covariation_partners$Protein_1,covariation_partners$Protein_2))
# number of proteins with at least 1 covariation Partner
universe |> length()
# universe are all proteins with at least 1 covariation partner
universe_geneID = Uniprot_annot[Entry %in%universe ,.(Entry,GeneID)]

covariation_partners = merge(covariation_partners,
                             Uniprot_annot[,.(Entry,GeneID)],
                             by.x = 'Protein_1',
                             by.y = 'Entry') |> 
  merge(Uniprot_annot[,.(Entry,GeneID)],
        by.x = 'Protein_2',
        by.y = 'Entry') 
covariation_partners = covariation_partners[GeneID.x != ''
][GeneID.y != ''   ]

enrichment_function <- function(Cluster_number,Universe){
  members  = clusterONE[Cluster==Cluster_number, Members]
  members = str_split(members,' ',simplify = T) |> unlist()
  members = universe_geneID[Entry %in% members,GeneID]
  if(length(members)>2){
    ego <- enrichGO(gene          =members ,
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
        term_OI[,Cluster:=Cluster_number]
        fwrite(term_OI, here::here('out','datasets',glue::glue('enrichment_terms_clusterONE.csv')), 
               append = T)
        
      }
      
    }
  }
}
library(future)
plan(sequential)
plan(multisession, workers = 20)

furrr::future_walk(.x = clusterONE$Cluster, ~enrichment_function(Cluster_number = .x,
                                                                    Universe = universe_geneID ))
cluster_enrichment = fread(here::here('out','datasets',glue::glue('enrichment_terms_clusterONE.csv')))
cluster_enrichment[, Cluster:=as.character(Cluster)]
Cluster_enrichment = merge(cluster_enrichment[order(p.adjust)][,.SD[head(1)], by = Cluster],
      clusterONE[to_remove==F], by ='Cluster',all.y = T)
fwrite(Cluster_enrichment, here::here('out','datasets','ClusterOne_clusters.gz'))
Cluster_enrichment = fread(here::here('out','datasets','ClusterOne_clusters.gz'))
is.na(Cluster_enrichment$ONTOLOGY) |> table()         

clusters_covar = data.table()
for(i in Cluster_enrichment$Cluster){
  print(i)
  cluster_members = Cluster_enrichment[Cluster == i,Members] |> 
    str_split(' ',simplify = T) |> as.vector()
  clusters_covar = rbind(clusters_covar,
                         bicor_covar[Protein_1 %in% cluster_members &
                                       Protein_2 %in% cluster_members 
                                     ][,.(avg_covar = mean(corr_bicor,na.rm = T),
                                          sd_covar = sd(corr_bicor,na.rm = T),
                                          N_protens = .N
                                          )][,Cluster := i])
}
Cluster_enrichment = merge(clusters_covar,Cluster_enrichment, by ='Cluster') 
  ggplot(Cluster_enrichment, aes(x = avg_covar, y = sd_covar, label =Description, size = Size ))+
  geom_point()+
  theme_bw()+
  ggrepel::geom_text_repel( data = Cluster_enrichment[avg_covar>0.6 & ! is.na(Description)],size = 3)+
    ggtitle('ClusterOne Clusters')
  ggsave(here::here('out','plots','ClusterONE_top_clusters.pdf'))
 
  
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
  
  
  bicorcor_OI_int_piv  = df |>
    as.data.frame() |> 
    tibble::rownames_to_column('Uniprot') |> 
    as.data.table() |> 
    melt(id.vars = 'Uniprot', variable.name = 'experiment',value.name = 'LogRatio')
  
  clusters_to_plot =  Cluster_enrichment[order(-avg_covar),]
  
  bicorcor_OI_int_piv[,condition:= str_extract(experiment,'_S[:print:]*') |> str_remove('^_')]
  bicorcor_OI_int_piv[,category:=dplyr::case_when(
    Uniprot %in% as.vector(str_split(clusters_to_plot[1,Members],' ',simplify = T)) ~clusters_to_plot[1,Description],
    Uniprot %in% as.vector(str_split(clusters_to_plot[2,Members],' ',simplify = T)) ~clusters_to_plot[2,Description],
    Uniprot %in% as.vector(str_split(clusters_to_plot[6,Members],' ',simplify = T)) ~clusters_to_plot[6,Description],
    Uniprot %in% as.vector(str_split(clusters_to_plot[5,Members],' ',simplify = T)) ~clusters_to_plot[5,Description],
    TRUE~'other'
  )]
  bicorcor_OI_int_piv =bicorcor_OI_int_piv[category!='other'
                                           ]
  colour_pallet = c('darkgreen','darkred','darkblue','darkorange')
  N_proteins = bicorcor_OI_int_piv[!is.na(LogRatio)][,.(N_quant = .N), by = experiment]
  experiments_to_plots = N_proteins[N_quant>quantile(N_proteins$N_quant)[2],experiment] |> 
    head(100)
  bicorcor_OI_int_piv = bicorcor_OI_int_piv[experiment %in% experiments_to_plots]
  ggplot(bicorcor_OI_int_piv,aes(x = experiment, y= LogRatio,
                                 colour= category,group=Uniprot))+
    # geom_line(data = bicorcor_OI_int_piv[category == 'other'],alpha = 0.05)+
    geom_line()+
    # geom_line(data = bicorcor_OI_int_piv[!(category %in% c(check_Protein_1,check_Protein_2))],alpha = 0.6)+
    # geom_line(data = bicorcor_OI_int_piv[(category %in% c(check_Protein_1,check_Protein_2))],alpha = 1)+
    theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),
          legend.position = 'bottom',
          panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    scale_colour_manual(values = colour_pallet)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    ggtitle('ClusterONE top covarying clusters')
  
  ggsave(here::here('out','plots','ClusterONE_top_covariation.pdf'))
  

  # between cluster covariation
  cluster_annotation = clusters_to_plot[!is.na(ID)][1:8,.(Members,Description)] |> 
    tidyr::separate_longer_delim(cols = Members,delim = ' ') |> 
    tibble::column_to_rownames('Members')
  
  cluster_members = bicor_covar[Protein_1 %in% rownames(cluster_annotation) &
                                  Protein_2 %in%  rownames(cluster_annotation)]

  cluster_members_rev = copy(cluster_members)
  cluster_members_rev[,`:=`(Protein_1= Protein_2,
                            Protein_2= Protein_1)]
  cluster_members = rbind(cluster_members_rev,
                          cluster_members) |> 
    dcast(Protein_1 ~ Protein_2, value.var = 'corr_bicor') 
  cluster_members = merge(cluster_members,
                          Uniprot_annot[,.(Entry, `Gene Names`)], by.x = 'Protein_1',
                          by.y = 'Entry')
  cluster_members[,Protein_1:=NULL]  
  setnames(cluster_members,'Gene Names','Protein_1')
  cluster_members = cluster_members |>  tibble::column_to_rownames('Protein_1') |> 
    as.matrix()
  diag(cluster_members) <- 1
  cluster_members[is.na(cluster_members)] <- rnorm(sum(is.na(cluster_members)),
                                                   0,0.1)
  pheatmap::pheatmap(cluster_members, 
                     # annotation_row =cluster_annotation,
                     annotation_col = cluster_annotation,
                     )
  
  # between cluster covariation
  cluster_annotation = clusters_to_plot[!is.na(ID)
                                        ][,Description:= paste(Description,Cluster,sep ='_')
                                          ][1:8,.(Members,Description)] |> 
    tidyr::separate_longer_delim(cols = Members,delim = ' ') |> 
    tibble::column_to_rownames('Members')
  
  cluster_members = bicor_covar[Protein_1 %in% rownames(cluster_annotation) &
                                  Protein_2 %in%  rownames(cluster_annotation)]
  
  cluster_members_rev = copy(cluster_members)
  cluster_members_rev[,`:=`(Protein_1= Protein_2,
                            Protein_2= Protein_1)]
  cluster_members = rbind(cluster_members_rev,
                          cluster_members) |> 
    dcast(Protein_1 ~ Protein_2, value.var = 'corr_bicor') 
  cluster_members = merge(cluster_members,
                          Uniprot_annot[,.(Entry, `Gene Names`)], by.x = 'Protein_1',
                          by.y = 'Entry')
  cluster_members[,Protein_1:=NULL]  
  setnames(cluster_members,'Gene Names','Protein_1')
  cluster_members = cluster_members |>  tibble::column_to_rownames('Protein_1') |> 
    as.matrix()
  diag(cluster_members) <- 1
  cluster_members[is.na(cluster_members)] <- rnorm(sum(is.na(cluster_members)),
                                                   0,0.1)
  pheatmap::pheatmap(cluster_members, 
                     # annotation_row =cluster_annotation,
                     annotation_col = cluster_annotation,
  )
  
  
  # finding similar complexes by function not by subunits or covariation
  clusters_to_group = cluster_enrichment[!(Cluster %in% clusters__to_remove)]
  list_clusters = list()
  for(i in unique(clusters_to_group$Cluster)){
    list_clusters[[i]] <- clusters_to_group[Cluster ==i,ID]
  }
  clusters_GO_similarity = calc_pairwise_overlaps(list_clusters) |> 
    as.data.table()
  clusters_GO_similarity_rev = copy(clusters_GO_similarity)
  clusters_GO_similarity_rev[,`:=`(name1= name2,
                            name2= name1)]
  clusters_GO_similarity = rbind(clusters_GO_similarity_rev,
                          clusters_GO_similarity) |> 
    dcast(name1 ~ name2, value.var = 'overlap') 
  clusters_GO_similarity = clusters_GO_similarity |>  tibble::column_to_rownames('name1') |> 
    as.matrix()
  diag(clusters_GO_similarity) <- 1
  pheatmap::pheatmap(clusters_GO_similarity,fontsize = 6
  )
  
  # between cluster covariation
  cluster_annotation = Cluster_enrichment[Cluster %in% c('232','242','305','135','226')
  ][,Description:= paste(Description,Cluster,sep ='_')
  ][,.(Members,Description)] |> 
    tidyr::separate_longer_delim(cols = Members,delim = ' ') |> 
    tibble::column_to_rownames('Members')
  
  cluster_members = bicor_covar[Protein_1 %in% rownames(cluster_annotation) &
                                  Protein_2 %in%  rownames(cluster_annotation)]
  
  cluster_members_rev = copy(cluster_members)
  cluster_members_rev[,`:=`(Protein_1= Protein_2,
                            Protein_2= Protein_1)]
  cluster_members = rbind(cluster_members_rev,
                          cluster_members) |> 
    dcast(Protein_1 ~ Protein_2, value.var = 'corr_bicor') 
  cluster_members = merge(cluster_members,
                          Uniprot_annot[,.(Entry, `Gene Names`)], by.x = 'Protein_1',
                          by.y = 'Entry')
  cluster_members[,Protein_1:=NULL]  
  setnames(cluster_members,'Gene Names','Protein_1')
  cluster_members = cluster_members |>  tibble::column_to_rownames('Protein_1') |> 
    as.matrix()
  diag(cluster_members) <- 1
  cluster_members[is.na(cluster_members)] <- rnorm(sum(is.na(cluster_members)),
                                                   0,0.1)
  pheatmap::pheatmap(cluster_members, 
                     # annotation_row =cluster_annotation,
                     annotation_col = cluster_annotation,
  )
  Cluster_enrichment_members = Cluster_enrichment |> 
    tidyr::separate_longer_delim(cols = 'Members',delim = ' ') |> 
    as.data.table()
  
  cluster_members_bicor =  bicor_covar[Protein_1 %in% Cluster_enrichment_members$Members & 
                Protein_2 %in% Cluster_enrichment_members$Members]
  cluster_members_bicor_rev = copy(cluster_members_bicor)
  cluster_members_bicor_rev[,`:=`(Protein_1 = Protein_2,
         Protein_2 = Protein_1)]
  cluster_members_bicor = rbind(cluster_members_bicor,
                                cluster_members_bicor_rev)
  cluster_members_bicor = merge(cluster_members_bicor, Cluster_enrichment_members[,.(Cluster,Members)], allow.cartesian=TRUE,
                                by.x = 'Protein_1', by.y = 'Members') |> 
    merge(Cluster_enrichment_members[,.(Cluster,Members)], allow.cartesian=TRUE,by.x = 'Protein_2', by.y = 'Members')
  avg_cluster_covar  = cluster_members_bicor[,.(cluster_pair_covar = mean(corr_bicor,na.rm = T)), 
                                             by =.(Cluster.x,Cluster.y)
                                             ]
  matrix_clusters = copy(avg_cluster_covar)
  matrix_clusters[,`:=`(Cluster.x = as.character(Cluster.x),
                        Cluster.y = as.character(Cluster.y))]
  matrix_clusters = matrix_clusters |> dcast(Cluster.x ~ Cluster.y, value.var = 'cluster_pair_covar') |> 
    tibble::column_to_rownames('Cluster.x')
  matrix_clusters = matrix_clusters[colnames(matrix_clusters),colnames(matrix_clusters)]
  diag(matrix_clusters) <- 1
  # matrix_clusters[is.na(matrix_clusters)] <- rnorm(sum(is.na(matrix_clusters)),0,0.1)
  pheatmap::pheatmap(matrix_clusters, cluster_rows = T, cluster_cols = T)
  
  to_join = clusters_overlap[,.(name1,name2,overlap)]
  setnames(to_join,c('Cluster.x','Cluster.y','overlap'))
  to_join[,`:=`(Cluster.x = as.integer(Cluster.x),
                Cluster.y = as.integer(Cluster.y))]
  avg_cluster_covar = merge(avg_cluster_covar,to_join, by = c('Cluster.x','Cluster.y'))
  
  avg_cluster_covar = avg_cluster_covar[Cluster.x>Cluster.y]  
  avg_cluster_covar |> 
    ggplot(aes(x = cluster_pair_covar , y = overlap))+
    geom_point(alpha = 0.5)+
    theme_bw()
  
  
  clusters_to_plot =  avg_cluster_covar[overlap==0][order(-cluster_pair_covar)
                                                    ][1,.(Cluster.x, Cluster.y)] |> 
    unlist() |> as.character()
  
  bicorcor_OI_int_piv  = df |>
    as.data.frame() |> 
    tibble::rownames_to_column('Uniprot') |> 
    as.data.table() |> 
    melt(id.vars = 'Uniprot', variable.name = 'experiment',value.name = 'LogRatio')
  
  
  bicorcor_OI_int_piv[,category:=NULL]
  bicorcor_OI_int_piv[,category:=dplyr::case_when(
  Uniprot %in% as.vector(str_split(Cluster_enrichment[Cluster == clusters_to_plot[1] ,Members],' ',simplify = T)) ~Cluster_enrichment[Cluster == clusters_to_plot[1],Description],
    Uniprot %in% as.vector(str_split(Cluster_enrichment[Cluster == clusters_to_plot[2] ,Members],' ',simplify = T)) ~Cluster_enrichment[Cluster == clusters_to_plot[2],Description],
    TRUE~'other'
  )]
  bicorcor_OI_int_piv =bicorcor_OI_int_piv[category!='other'
  ]
  colour_pallet = c('darkgreen','darkorange')
  N_proteins = bicorcor_OI_int_piv[!is.na(LogRatio)][,.(N_quant = .N), by = experiment]
  experiments_to_plots = N_proteins[N_quant>=quantile(N_proteins$N_quant)[2],experiment] |>
    head(50)
  bicorcor_OI_int_piv = bicorcor_OI_int_piv[experiment %in% experiments_to_plots]
  ggplot(bicorcor_OI_int_piv,aes(x = experiment, y= LogRatio,
                                 colour= category,group=Uniprot))+
    # geom_line(data = bicorcor_OI_int_piv[category == 'other'],alpha = 0.05)+
    # geom_line(alpha = 0.5)+
    geom_line(data = bicorcor_OI_int_piv[str_detect(category,'regulation of gene expression')],alpha = 0.6)+
    geom_line(data = bicorcor_OI_int_piv[str_detect(category,'regulation of gene expression',negate = T)],alpha = 0.6)+
    theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),
          legend.position = 'bottom',
          panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    scale_colour_manual(values = colour_pallet)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    ggtitle('ClusterONE top covarying clusters',
            subtitle = 'have no shared members')
  
  ggsave(here::here('out','plots','ClusterONE_covarying_modules.pdf'))
``  
  
  