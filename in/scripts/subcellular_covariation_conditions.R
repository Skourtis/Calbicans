# for each subcelular compartment
# plot pairwise covariation between the proteins of the compartmetn
# and between compartments

DT_normalised = fread(here::here('out',
                                 'datasets',
                                 'covariation_stressed_treatments_normalised.gz'))
# Uniprot mapping
Uniprot_annot = fread(here::here('in','datasets','uniprotkb_taxonomy_id_237561_2025_07_08.tsv.gz'))
Uniprot_annot[,CGD:=str_remove(CGD,';[:print:]*$')]
uniprot_subcell = Uniprot_annot[,.(Entry,`Subcellular location [CC]`)]
uniprot_subcell[,localisation := str_remove_all(`Subcellular location [CC]`,'SUBCELLULAR LOCATION\\:')]
uniprot_subcell[,`Subcellular location [CC]`:=NULL]
uniprot_subcell[,localisation:= str_split(localisation,pattern = '\\}(.|;) ')]
uniprot_subcell = tidyr::unnest(uniprot_subcell,localisation,) |> as.data.table() |> unique()
uniprot_subcell[,localisation:= str_remove_all(localisation,' \\{[:print:]*$') |> trimws()]
uniprot_subcell_n = uniprot_subcell[,.(N_prots = .N), by = localisation]
uniprot_subcell[,N_prots := .N, by = localisation]
uniprot_subcell[localisation!='',final_localisation:=fifelse(N_prots>100,localisation,NA_character_)]
# uniprot_subcell[final_localisation %in% c('Cytoplasm','Secreted'),final_localisation:=NA_character_]


compartment_covariation = data.table()
for(compartment in na.omit(unique(uniprot_subcell$final_localisation))){
  print(compartment)
  compartment_proteins = uniprot_subcell[final_localisation  == compartment,Entry]
  compartment_covariation = rbind(compartment_covariation,
                                  DT_normalised[Protein_1 %in% compartment_proteins & 
                                                  Protein_2 %in% compartment_proteins
                                                ][,`:=`(final_localisation=compartment ,
                                                                 bands = NULL,
                                                                 N_pairs_bands = NULL)])
  
}
compartment_covariation |> 
  ggplot(aes(y = Condition, x = zscore))+
  geom_boxplot()+
  theme_bw()+
  facet_wrap('final_localisation')

# possible PPI complexes
consistent_interactors = DT_normalised[avg_cor>0.5,.(Protein_1,Protein_2,avg_cor,sd_cor)] |> 
  unique()
  ggplot(consistent_interactors,aes(x = avg_cor,y = sd_cor,
                                    label = paste(Protein_1,Protein_2,sep ='_')))+
  geom_point(alpha = 0.5)+
  theme_bw()+
    ggrepel::geom_text_repel(data = consistent_interactors[sd_cor<0.1 & avg_cor>0.8])
  
  
  
clusterone_modules =  fread(here::here('out','datasets','ClusterOne_clusters.gz'))
clusterone_modules[,Cluster := as.factor(Cluster)]
module_covariation = data.table()
for(module in clusterone_modules$Cluster){
  print(module)
  module_proteins = clusterone_modules[Cluster  == module,Members] |> 
    str_split(' ',simplify = T) |> as.vector()
  
  module_covariation = rbind(module_covariation,
                                  DT_normalised[Protein_1 %in% module_proteins & 
                                                  Protein_2 %in% module_proteins
                                  ][,`:=`(Cluster=module ,
                                          bands = NULL,
                                          N_pairs_bands = NULL)])
  
}
module_covariation_per_condition = module_covariation[,.(module_covariation = median(bicor,na.rm =T),
                                                         sd_condition = median(bicor,na.rm =T),
                                                         N_pairs = .N,
                                                         members = uniqueN(c(Protein_1,Protein_2))),by = .(Cluster,Condition)]

module_covariation_across_condition = module_covariation_per_condition[,.(avg_covariation = median(module_covariation),
                                                                          spreadness_within_conditions = median(sd_condition),
                                                                          sd_across_conditions = sd(module_covariation),
                                                                          N_pairs = mean(N_pairs),
                                                                          members = mean(members)),by = Cluster]
module_covariation_across_condition |> ggplot(aes(x =avg_covariation,
                                                  y = (sd_across_conditions),
                                                  size = members,
                                                  colour = (spreadness_within_conditions),
                                                  label = Cluster))+
  geom_point()+
  theme_bw()+
  scico::scale_color_scico()+
  ggrepel::geom_text_repel(size = 4)

examples = c('53','237','130','293')
examples_to_plot = module_covariation[ Cluster %in% examples] |> 
  merge(clusterone_modules[,.(Cluster, Description)], by = 'Cluster')
examples_to_plot[,Cluster_Description:= paste0('Cluster ',Cluster,'\n',Description)]
  ggplot(examples_to_plot,aes(x = Condition, y = bicor))+
  geom_jitter(height   = 0, width = 0.2, alpha = 0.1)+
  geom_boxplot(outliers = F, fill =NA)+
  facet_wrap('Cluster_Description')+
  theme_bw()+
    ggtitle('Examples of modules with their protein pairwise bicor members')
