# differential expression of modules
library('EnrichIntersect')
Cluster_one_modules = fread(here::here('out','datasets','ClusterOne_clusters.gz'))
modules_members = Cluster_one_modules[Size>6,.(Cluster,Members)] |> 
  tidyr::separate_longer_delim(cols = 'Members',delim = ' ') |> 
  as.data.table()
modules_members[,.N, by = Members][,N] |> 
  hist(main = 'most proteins in single module')

# Uniprot mapping
Uniprot_annot = fread(here::here('in','datasets','uniprotkb_taxonomy_id_237561_2025_07_08.tsv.gz'))
Uniprot_annot[,CGD:=str_remove(CGD,';[:print:]*$')]
# Uniprot_annot =Uniprot_annot[GeneID != '']

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
subset_stress = merge(subset_stress,Uniprot_annot[,.(Entry,CGD)],by.x = 'Genes',by.y = 'CGD')
subset_stress[,Genes:=NULL]
subset_stress <- subset_stress |> tibble::column_to_rownames('Entry')
subset_stress |> dim()

# Transpose and median normalise the rows to avoid spurious correlations
tmp_medians <- apply( subset_stress, 1, median, na.rm = TRUE )  
df <- sweep( subset_stress, 1, tmp_medians, FUN = "-" )


### complexes which are downregulated in specific tissues #### 
stressed_matrix_norm = copy(df) |> t()
prots_NA= colMeans(is.na(stressed_matrix_norm)) 
hist(prots_NA)
stressed_matrix_norm = stressed_matrix_norm[,prots_NA<0.8]
stressed_matrix_norm |> dim()
stressed_matrix_norm[,1:20] |> boxplot()
stressed_matrix_norm |> hist(breaks = 100)
# imputing the missing values of each protein to -4
stressed_matrix_norm[is.na(stressed_matrix_norm)] <- rnorm(sum(is.na(stressed_matrix_norm)),
                                                                       mean = -2,
                                                                       sd = 0.3)

tmp_medians <- apply( stressed_matrix_norm , 2, median, na.rm = TRUE )  
stressed_matrix_norm <- sweep( stressed_matrix_norm , 2,tmp_medians, FUN = "-" )

stressed_matrix_norm |> hist(breaks = 100)


complex_candida = stressed_matrix_norm |> as.data.frame() |> 
  tibble::rownames_to_column('Project_Identifier') |>
  as.data.table() |> 
  melt(id.vars = 'Project_Identifier', value.name = 'zscore_abundance',variable.name = 'Protein')
complex_candida |> nrow()
complex_candida = merge(complex_candida,  modules_members,by.x = 'Protein',by.y = 'Members')
complex_candida |> nrow()
complex_candida[,Condition:=str_extract(Project_Identifier,'_S[:print:]*$') |> str_remove('^_')]
complex_candida[,.(N_samples = uniqueN(Project_Identifier)),by =Condition ]

sample_mapping = complex_candida[,.(Project_Identifier,Condition)] |> 
  unique()

complex_candida_avg = complex_candida[,.(SD_module = sd(zscore_abundance),
                                              mean = mean(zscore_abundance),
                                              N_subunits = .N),by = .(Condition,Cluster)] 

complex_candida_enrichment = complex_candida[,.(SD_module = sd(zscore_abundance),
                                                         mean = mean(zscore_abundance),
                                                         N_subunits = .N),by = .(Condition,Project_Identifier , Cluster)] 

complex_candida_enrichment  = complex_candida_enrichment |>
  dcast(Project_Identifier ~ Cluster ,value.var = 'mean') |> 
  tibble::column_to_rownames('Project_Identifier') |> 
  as.matrix()


custom.set = data.table(drug = sample_mapping$Project_Identifier,
                        group = sample_mapping$Condition)
custom.set = custom.set[drug %in% rownames(complex_candida_enrichment)]

set.seed(123)
enrich <- enrichment(complex_candida_enrichment, custom.set, permute.n = 100)
values_complexes = enrich$pvalue |> as.data.frame() |> 
  tibble::rownames_to_column('Condition') |> 
  as.data.table() |> 
  melt(id.vars = 'Condition',value.name = 'pvalue', variable.name = 'Cluster')
values_complexes[,padj:= p.adjust(pvalue,method = 'BH')]
NES_complexes = enrich$S |> as.data.frame() |> 
  tibble::rownames_to_column('Condition') |> 
  as.data.table() |> 
  melt(id.vars = 'Condition',value.name = 'NES', variable.name = 'Cluster')
enrich_complexes = merge(NES_complexes,values_complexes, by = c('Condition','Cluster'))
enrich_complexes[,pvalue_plot := pvalue+rnorm(1,0.03,0.0065), by = .(Cluster,Condition) ]
# order_types = enrich_complexes[Cluster %in% c('huMAP3_00154.1')][order(-NES),Condition]
# enrich_complexes[,Condition := factor(Condition,levels =order_types)]
enrich_complexes |> 
  ggplot(aes(x = NES, y = -log10(pvalue_plot), label = Cluster))+
  geom_point(aes(alpha = padj<0.01))+
  theme_bw()+
  geom_point(data = enrich_complexes[Cluster %in% c('192',
                                                    '2',
                                                    '110')],
             aes(colour = Cluster))+
  ggrepel::geom_text_repel(data = enrich_complexes[Cluster %in% c('192',
                                                           '2',
                                                           '110')],
                           aes(colour = Cluster))+
  facet_wrap('Condition')+
  theme(legend.position = 'bottom')+
  labs(x = 'Normalised Enrichment Score', y = '-log10(pvalue)')
ggsave(here::here('out','plots','NES_Cluster_conditions.pdf'),width = 25,height = 15)


diff_complexes =   complex_candida_avg[mean>1 &SD_module<1.2][,.N,by = Cluster][order(N)] |> tail(3)
diff_complexes = diff_complexes[,Cluster]

original_data = subset_stress |> t() |> as.data.frame() |> 
  tibble::rownames_to_column('Project_Identifier') |>
  as.data.table() |> 
  melt(id.vars = 'Project_Identifier', value.name = 'imputation',variable.name = 'Protein')
plot_diff_complexes =  stressed_matrix_norm |> as.data.frame() |> 
  tibble::rownames_to_column('Project_Identifier') |>
  as.data.table() |> 
  melt(id.vars = 'Project_Identifier', value.name = 'abundance',variable.name = 'Protein')
plot_diff_complexes = merge(plot_diff_complexes,original_data, by = c('Project_Identifier','Protein'),
                            all = T)
plot_diff_complexes[,imputation:= fifelse(is.na(imputation),'imputted','not_imputted')]
plot_diff_complexes |> ggplot(aes(x = abundance, fill = imputation ))+
  geom_density(alpha = 0.5)

plot_diff_complexes = merge(plot_diff_complexes,  modules_members,
                            by.x = 'Protein', by.y = 'Members')
plot_diff_complexes[,Condition:=str_extract(Project_Identifier,'_S[:print:]*$') |> str_remove('^_')]
# plot_diff_complexes = merge(plot_diff_complexes, Gene_names, by = 'Protein')

# order_samples = complex_candida[Cluster == tail(diff_complexes,1)][order(mean),Condition]

# plot_diff_complexes[,Condition:=factor(Condition,levels = order_samples)]

plot_diff_complexes[1:10000,] |> 
  ggplot(aes(y = Condition, x = abundance ))+
  geom_boxplot()+
  geom_jitter( alpha = 0.1,aes(colour = imputation))+
  facet_wrap('Cluster')+theme_bw()+
  labs(x = 'Protein Relative Abundance', colour = 'Imputation', y = 'Cancer Lineage')
  ggrepel::geom_text_repel(data = plot_diff_complexes[order(abundance),head(.SD,1),by = .(ID,Cluster)],
                           max.overlaps = 40,aes(y = Condition, x = abundance,label = ID))

ggsave(here::here('out','plots','diff_complex_abundance.pdf'),width=9, height = 9)
fwrite(plot_diff_complexes,here::here('out','plots','diff_complex_abundance.gz'))

order_to_plot = plot_diff_complexes[Cluster == '192',
                                    ][,.(avg_expression = mean(abundance)),by = Condition
                                      ][order(avg_expression),Condition ]
example_to_plot = plot_diff_complexes[Cluster %in% c('192','2','110'),
                                      ][,Condition := factor(Condition,
                                                             levels = order_to_plot)]
example_to_plot  = merge(example_to_plot,
                         Cluster_one_modules[,.(Cluster,Description)],
                         by = 'Cluster')
example_to_plot[,Cluster_Description:= paste0('Cluster ',Cluster,'\n',Description)]
example_to_plot|> 
  ggplot(aes(y = Condition, x = abundance ))+
  geom_jitter( alpha = 0.1,aes(colour = imputation))+
  geom_boxplot(outliers = F, fill =NA)+
  facet_wrap('Cluster_Description')+theme_bw()+
  labs(x = 'Protein Relative Abundance', colour = 'Imputation', y = 'Cancer Lineage')

