# Load necessary libraries
library(data.table);library(ggplot2); library(gridExtra); library(umap);library(stringr)
set.seed(123)
# Uniprot mapping
Uniprot_annot = fread(here::here('in','datasets','uniprotkb_taxonomy_id_237561_2025_07_07.tsv.gz'))
Uniprot_annot[,CGD:=str_remove(CGD,';[:print:]*$')]
#### Prepare co-regulation scores ####
DT = fread(here::here('out','datasets','whole_streesed_covariations.gz'))
# Turn co-regulation score back into a "distance" metric and log2-transform for better tSNE performance
DT[, `:=`  (coreg_distance = (1 - (corr_bicor+1)/2))]

# Turn the melted pairwise table back into a dist object
DTm <- dcast( data = rbind( DT[, .(Protein_1, Protein_2, coreg_distance)],                             # These steps create a "redundant" table...
                            DT[, .(Protein_1 = Protein_2, Protein_2 = Protein_1, coreg_distance)]),    # ... containing both A <-> B and B <-> pairs
              Protein_1 ~ Protein_2 , value.var = "coreg_distance")                      # And this casts them into a matrix-shaped data.table
DTm <- as.data.frame(DTm)                 # Turn into data.frame
rownames(DTm) <- DTm$Protein_1            # Create rownames
DTm$Protein_1 <- NULL                     # Drop original name column
DTm <- as.matrix(DTm)

# adding some noise to RF scores
DTm_noise = DTm 
# + matrix(rnorm(ncol(DTm)^2,0,0.0015),nrow = ncol(DTm),ncol = ncol(DTm))
diag(DTm_noise) <- 0

# imputting missing pairwwise with low values anre makign the imputation symmetrical
DTm_noise[is.na(DTm_noise)] <- rbeta(is.na(DTm_noise) |> sum(),8.5,2.2)
DT_m_sym = Matrix::forceSymmetric(DTm_noise |> as.matrix())
DT_m_sym = as.matrix(DT_m_sym)

umap_DT_coreg_umap= umap::umap(DT_m_sym, input = 'dist')
#### Create umap map ####
uMap_layout <- data.table( umap_DT_coreg_umap$layout)
uMap_layout[, ID := rownames(umap_DT_coreg_umap$data )]
names(uMap_layout) <- c("UMAP1", "UMAP2", "ID")
fwrite(uMap_layout[,.(Entry,UMAP1,UMAP2)],
       here::here('out','datasets','UMAP_covariation.gz'))
# uMap_layout =umap_DT_coreg_uwot |> as.data.frame() |> tibble::rownames_to_column('ID') |> as.data.table()
# names(uMap_layout) <- c( "ID","UMAP1", "UMAP2")
pumap = ggplot(uMap_layout, aes(x = UMAP1, y = UMAP2))+
  geom_point(shape = 16, size=0.9, alpha=0.8)+
  # geom_point(data = uMap_layout[between(UMAP1,2,3) & between(UMAP2,-2,-0.5)], colour ='red')+
  theme(panel.background=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(),
        panel.border=element_rect(fill=NA, colour="black", size=0.95), panel.grid.major=element_blank(),
        axis.title=element_blank(), plot.margin = unit( c(0.5, 0.5, 0.5, 0.5), "mm"))

ggsave(here::here('out','plots','bicor_umap_no_annot.pdf'),pumap)
# Load and pre-process the yeast GO slim, which was downloaded from:
uniprot_subcell = Uniprot_annot[,.(Entry,`Subcellular location [CC]`)]
uniprot_subcell[,localisation := str_remove_all(`Subcellular location [CC]`,'SUBCELLULAR LOCATION\\:')]
uniprot_subcell[,`Subcellular location [CC]`:=NULL]
uniprot_subcell[,localisation:= str_split(localisation,pattern = '\\}(.|;) ')]
uniprot_subcell = tidyr::unnest(uniprot_subcell,localisation,) |> as.data.table() |> unique()
uniprot_subcell[,localisation:= str_remove_all(localisation,' \\{[:print:]*$') |> trimws()]
uniprot_subcell_n = uniprot_subcell[,.(N_prots = .N), by = localisation]
uniprot_subcell[,N_prots := .N, by = localisation]
uniprot_subcell[localisation!='',final_localisation:=fifelse(N_prots>100,localisation,NA_character_)]
uniprot_subcell[final_localisation %in% c('Cytoplasm','Secreted'),final_localisation:=NA_character_]
# uniprot_subcell[,final_localisation:= fcase(
#   str_detect(localisation,'^Multi-pass'),'Multi-pass membrane protein',
#   str_detect(localisation,'ucleol'),'Nucleolus',
#   str_detect(localisation,'ucleus'),'Nucleus',
#   str_detect(localisation,'itochon'),'Mitochondria',
#   str_detect(localisation,'eticulu'),'ER',
#   str_detect(localisation,'olgi'),'Golgi',
#   str_detect(localisation,'ecreted'),'Secreted',
#   str_detect(localisation,'ytoplas'),'Cytoplasm'
# ) ]

# uniprot_subcell = uniprot_subcell[!is.na(final_localisation)]
uniprot_subcell[,N_terms_prot := .N, by = Entry]
uniprot_subcell[,N_terms_comp := .N, by = Entry]
uniprot_subcell_N = uniprot_subcell[,.(N_prots = .N), by = localisation]
# Set plotting colours
col_values = c( "Cytoplasm" = "#b5c900",
                "Nucleus" = "#0044e3", 
                "Nucleolus" = "#009fe3",
                "Secreted" = "#00e3cc",
                "Multi-pass membrane protein" = "#9700e3",
                "Mitochondria" = "#ff0080",
                "ER" = "#e3b200",
                "Golgi" = "#ffee00")
uMap_layout[,Entry := str_remove(ID,'\\.[:print:]*$')]
uniprot_subcell_plot = merge(uMap_layout,uniprot_subcell, by = 'Entry', all.x = T )
scico::scico_palette_names(categorical = T)
# Plot the subcellular locations on the maps
ggplot(uniprot_subcell_plot, aes(UMAP1, UMAP2, colour = final_localisation))+
  # facet_wrap('final_localisation', scale = "free")+
  geom_point(data = uniprot_subcell_plot[is.na(final_localisation)],  size = 0.5, alpha = 0.6, shape = 16, colour = "grey70")+
  geom_point(data = uniprot_subcell_plot[!is.na(final_localisation)],  size = 1, alpha = 0.8)+
  # # geom_point(size = 0.8, alpha = 0.9, shape = 16, data = DT[ SimpleGO == "Cytoplasm" ])+
  
  # geom_point(size = 1, alpha = 0.9, shape = 16, data = uniprot_subcell_plot[ final_localisation == "Nucleus" ])+
  # geom_point(size = 1, alpha = 0.9, shape = 16, data = uniprot_subcell_plot[ final_localisation == "Cytoplasm" ])+
  # geom_point(size = 1, alpha = 0.9, shape = 16, data = uniprot_subcell_plot[ final_localisation == "Nucleolus" ])+
  # geom_point(size = 1, alpha = 0.9, shape = 16, data = uniprot_subcell_plot[ final_localisation == "Secreted" ])+
  # geom_point(size = 1, alpha = 0.9, shape = 16, data = uniprot_subcell_plot[ final_localisation == "ER" ])+
  # geom_point(size = 1, alpha = 0.9, shape = 16, data = uniprot_subcell_plot[ final_localisation == "Golgi" ])+
  # geom_point(size = 1, alpha = 0.9, shape = 16, data = uniprot_subcell_plot[ final_localisation == "Multi-pass membrane protein" ])+
  # geom_point(size = 1, alpha = 0.9, shape = 16, data = uniprot_subcell_plot[ final_localisation == "Mitochondria" ])+
  # scale_colour_manual(values = col_values, name = "Subcellular location")+
  theme_bw()+ theme(legend.position="bottom")+
  scico::scale_color_scico_d(palette = 'batlow')+
  facet_wrap('final_localisation')
ggsave( here::here('out','plots', 'ProHD_RF_uniprot_localisation.pdf'))

