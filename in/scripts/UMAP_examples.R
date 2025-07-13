# Load necessary libraries
# install.packages('')
library(data.table);library(ggplot2); library(gridExtra); library(umap); library(stringr)
Uniprot_mapping =fread(here::here('in','datasets','uniprotkb_taxonomy_id_237561_2025_07_07.tsv.gz'))
Uniprot_mapping[,ID:=str_remove(`Gene Names`,' [:print:]*$')]
Uniprot_mapping[,CGD:=str_remove(CGD,';[:print:]*$')]
set.seed(1234)
# Load the coregulation scores (treeClust + TOM)
DT = fread(here::here('out','datasets','whole_streesed_covariations.gz'))
DT = DT[corr_bicor>=0.507,.(Protein_1,Protein_2,corr_bicor)]
#loading UMAP from previous script 
uMap_layout = fread(here::here('out','datasets','UMAP_covariation.gz'))
uMap_layout = merge(uMap_layout,
      Uniprot_mapping[,.(Entry,ID)],by = 'Entry')
# selected proteins to check
to_check = examples_to_show[`Protein existence`=='Predicted' &
                              Annotation==1][(order(min_padj))][25:28,Entry]
to_check =  c('A0A1D8PCX5','A0A1D8PD99','A0A1D8PTA4','Q5ALU1','Q59Z69')[1:4]
uMap_layout[,annot:= fcase(
  Entry ==to_check[1] ,to_check[1],
  Entry ==to_check[2] ,to_check[2],
  Entry ==to_check[3] ,to_check[3],
  Entry ==to_check[4] ,to_check[4],
  Entry %in% unlist(DT[Protein_1 ==to_check[1] |Protein_2 ==to_check[1] ][,.(Protein_1,Protein_2)]),as.character(glue::glue('Int_{to_check[1]}')),
  Entry %in% unlist(DT[Protein_1 ==to_check[2] |Protein_2 ==to_check[2] ][,.(Protein_1,Protein_2)]),as.character(glue::glue('Int_{to_check[2]}')),
  Entry %in% unlist(DT[Protein_1 ==to_check[3] |Protein_2 ==to_check[3] ][,.(Protein_1,Protein_2)]),as.character(glue::glue('Int_{to_check[3]}')),
  Entry %in% unlist(DT[Protein_1 ==to_check[4] |Protein_2 ==to_check[4] ][,.(Protein_1,Protein_2)]),as.character(glue::glue('Int_{to_check[4]}'))
  
)]
colour_annot_names = c(to_check[1],
                       to_check[2],
                       to_check[3],
                       to_check[4],
                       glue::glue('Int_{to_check[1]}'),
                       glue::glue('Int_{to_check[2]}') ,
                       glue::glue('Int_{to_check[3]}') ,
                       glue::glue('Int_{to_check[4]}'))
colour_annot = c('#228B22' ,
                 '#00008B',
                 '#CC5500',
                 'darkred',
                 '#32CD32',
                 '#89CFF0' ,
                  '#E3963E',
                 'pink')
names(colour_annot) = colour_annot_names
xcoord = uMap_layout[Entry %in% to_check,UMAP1]
ycoord = uMap_layout[Entry %in% to_check,UMAP2]
ggplot(uMap_layout, aes(x = UMAP1, y = UMAP2, colour =annot, label = ID ))+
  annotate('rect', xmin = xcoord[1]-0.3,xmax =  xcoord[1]+0.3, ymin = ycoord[1]-0.3, ymax = ycoord[1]+0.3, linewidth = 0.5, fill = 'grey85', colour= 'grey70')+
  annotate('rect', xmin = xcoord[2]-0.3,xmax =  xcoord[2]+0.3, ymin = ycoord[2]-0.3, ymax = ycoord[2]+0.3, linewidth = 0.5, fill = 'grey85', colour= 'grey70')+
  annotate('rect', xmin = xcoord[3]-0.3,xmax =  xcoord[3]+0.3, ymin = ycoord[3]-0.3, ymax = ycoord[3]+0.3, linewidth = 0.5, fill = 'grey85', colour= 'grey70')+
  annotate('rect', xmin = xcoord[4]-0.3,xmax =  xcoord[4]+0.3, ymin = ycoord[4]-0.3, ymax = ycoord[4]+0.3, linewidth = 0.5, fill = 'grey85', colour= 'grey70')+
  geom_point(shape = 16, size=0.9, alpha=0.2)+
  
  geom_point(data = uMap_layout[!is.na(annot)][str_detect(annot,'^Int')], shape = 16, size=1.2, alpha=0.8)+
  geom_point(data = uMap_layout[!is.na(annot)][str_detect(annot,'^Int',negate = T)], shape = 16, size=1.5, alpha=1)+
  theme(panel.background=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), legend.position = 'none',
        panel.border=element_rect(fill=NA, colour="black", size=0.95), panel.grid.major=element_blank(),
        axis.title=element_blank(), plot.margin = unit( c(0.5, 0.5, 0.5, 0.5), "mm"))+
  scale_color_manual(values = colour_annot)+
  ggrepel::geom_text_repel(data = uMap_layout[!is.na(annot)][str_detect(annot,'^Int',negate = T)])
ggsave(  here::here('out','plots','umap_poor_characterised_prots_interactors.pdf'),p1)

# statistical Go enrichments terms for significant functions
enrichments_bicor = fread(here::here('out','datasets',glue::glue('enrichment_terms_bicor.csv')))

# for each protein of interest zoom in
term_OI = enrichments_bicor[Protein == to_check[1]][Count>1][order(p.adjust)]
Interactors = term_OI[2,"geneID"] |> unlist() |> strsplit('/') |> unlist()

umapx = uMap_layout[Entry %in% to_check[1],UMAP1]
umapy = uMap_layout[Entry %in% to_check[1],UMAP2]
uMap_layout[,term:= fifelse(str_detect(ID, paste(Interactors,collapse = '|'),
                                       ),term_OI[2,Description],NA_character_)]
ggplot(uMap_layout, aes(x = UMAP1, y = UMAP2,  label = ID))+
  geom_point(data = uMap_layout[is.na(term ) ], shape = 16, size=3, alpha=0.7, colour = 'grey70')+
  geom_point(data = uMap_layout[!is.na(term ) ],shape = 16, size=3, alpha=0.7, colour =colour_annot[5] )+
  
  geom_point(data = uMap_layout[Entry == to_check[1]],shape = 16, size=3, alpha=1, colour = colour_annot[1])+
  # geom_point(data = uMap_layout[!is.na(annot)][str_detect(annot,'^Int')], shape = 16, size=1, alpha=0.9)+
  # geom_point(data = uMap_layout[!is.na(annot)][str_detect(annot,'^Int',negate = T)], shape = 16, size=1, alpha=0.9)+
  theme(panel.background=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), legend.position = 'none',
        panel.border=element_rect(fill=NA, colour="black", size=0.95), panel.grid.major=element_blank(),
        axis.title=element_blank(), plot.margin = unit( c(0.5, 0.5, 0.5, 0.5), "mm"))+
  ggrepel::geom_text_repel(data = uMap_layout[!is.na(term ) ], colour  =colour_annot[1])+
  ggrepel::geom_label_repel(data = uMap_layout[Entry == to_check[1]], fill =colour_annot[1], colour = 'white' )+
  # scale_color_manual(values = c('grey50',colour_annot[4]))+
  # lims(y = c())
  scale_x_continuous(limits = c(umapx-0.30,umapx+0.3),expand = c(0,0))+
  scale_y_continuous(limits = c(umapy-0.3,umapy+0.3),expand = c(0,0))+
  ggtitle(glue::glue('{to_check[1]} {term_OI[2,Description]}'),
          subtitle = glue::glue('p.adj {format(term_OI[2,p.adjust],nsmall = 3)}')) 
ggsave(here::here('out','plots','umap_closeup_1.pdf'))

term_OI = enrichments_bicor[Protein == to_check[2]][Count>1][order(p.adjust)]
Interactors = term_OI[2,"geneID"] |> unlist() |> strsplit('/') |> unlist()

umapx = uMap_layout[Entry %in% to_check[2],UMAP1]
umapy = uMap_layout[Entry %in% to_check[2],UMAP2]
uMap_layout[,term:= fifelse(str_detect(ID, paste(Interactors,collapse = '|'),
),term_OI[2,Description],NA_character_)]
ggplot(uMap_layout, aes(x = UMAP1, y = UMAP2,  label = ID))+
  geom_point(data = uMap_layout[is.na(term ) ], shape = 16, size=3, alpha=0.7, colour = 'grey70')+
  geom_point(data = uMap_layout[!is.na(term ) ],shape = 16, size=3, alpha=0.7, colour =colour_annot[6] )+
  
  geom_point(data = uMap_layout[Entry == to_check[2]],shape = 16, size=3, alpha=1, colour = colour_annot[2])+
  # geom_point(data = uMap_layout[!is.na(annot)][str_detect(annot,'^Int')], shape = 16, size=1, alpha=0.9)+
  # geom_point(data = uMap_layout[!is.na(annot)][str_detect(annot,'^Int',negate = T)], shape = 16, size=1, alpha=0.9)+
  theme(panel.background=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), legend.position = 'none',
        panel.border=element_rect(fill=NA, colour="black", size=0.95), panel.grid.major=element_blank(),
        axis.title=element_blank(), plot.margin = unit( c(0.5, 0.5, 0.5, 0.5), "mm"))+
  ggrepel::geom_text_repel(data = uMap_layout[!is.na(term ) ], colour  =colour_annot[6])+
  ggrepel::geom_label_repel(data = uMap_layout[Entry == to_check[2]], fill =colour_annot[2], colour = 'white' )+
  # scale_color_manual(values = c('grey50',colour_annot[4]))+
  # lims(y = c())
  scale_x_continuous(limits = c(umapx-0.3,umapx+0.3),expand = c(0,0))+
  scale_y_continuous(limits = c(umapy-0.3,umapy+0.3),expand = c(0,0))+
  ggtitle(glue::glue('{to_check[2]} {term_OI[2,Description]}'),
          subtitle = glue::glue('p.adj {format(term_OI[2,p.adjust],nsmall = 3)}')) 
ggsave(here::here('out','plots','umap_closeup_2.pdf'))

term_OI = enrichments_bicor[Protein == to_check[3]][Count>1][order(p.adjust)]
Interactors = term_OI[2,"geneID"] |> unlist() |> strsplit('/') |> unlist()

umapx = uMap_layout[Entry %in% to_check[3],UMAP1]
umapy = uMap_layout[Entry %in% to_check[3],UMAP2]
uMap_layout[,term:= fifelse(str_detect(ID, paste(Interactors,collapse = '|'),
),term_OI[2,Description],NA_character_)]
  
ggplot(uMap_layout, aes(x = UMAP1, y = UMAP2,  label = ID))+
  geom_point(data = uMap_layout[is.na(term ) ], shape = 16, size=3, alpha=0.7, colour = 'grey70')+
  geom_point(data = uMap_layout[!is.na(term ) ],shape = 16, size=3, alpha=0.7, colour =colour_annot[7] )+
  
  geom_point(data = uMap_layout[Entry == to_check[3]],shape = 16, size=3, alpha=1, colour = colour_annot[3])+
  # geom_point(data = uMap_layout[!is.na(annot)][str_detect(annot,'^Int')], shape = 16, size=1, alpha=0.9)+
  # geom_point(data = uMap_layout[!is.na(annot)][str_detect(annot,'^Int',negate = T)], shape = 16, size=1, alpha=0.9)+
  theme(panel.background=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), legend.position = 'none',
        panel.border=element_rect(fill=NA, colour="black", size=0.95), panel.grid.major=element_blank(),
        axis.title=element_blank(), plot.margin = unit( c(0.5, 0.5, 0.5, 0.5), "mm"))+
  ggrepel::geom_text_repel(data = uMap_layout[!is.na(term ) ], colour  =colour_annot[7])+
  ggrepel::geom_label_repel(data = uMap_layout[Entry == to_check[3]], fill =colour_annot[3], colour = 'white' )+
  # scale_color_manual(values = c('grey50',colour_annot[4]))+
  # lims(y = c())
  scale_x_continuous(limits = c(umapx-0.3,umapx+0.3),expand = c(0,0))+
  scale_y_continuous(limits = c(umapy-0.3,umapy+0.3),expand = c(0,0))+
  ggtitle(glue::glue('{to_check[3]} {term_OI[2,Description]}'),
          subtitle = glue::glue('p.adj {format(term_OI[2,p.adjust],nsmall = 3)}')) 
ggsave(here::here('out','plots','umap_closeup_3.pdf'))


term_OI = enrichments_bicor[Protein == to_check[4]][Count>1][order(p.adjust)]
Interactors = term_OI[2,"geneID"] |> unlist() |> strsplit('/') |> unlist()

umapx = uMap_layout[Entry %in% to_check[4],UMAP1]
umapy = uMap_layout[Entry %in% to_check[4],UMAP2]
uMap_layout[,term:= fifelse(str_detect(ID, paste(Interactors,collapse = '|'),
),term_OI[2,Description],NA_character_)]

ggplot(uMap_layout, aes(x = UMAP1, y = UMAP2,  label = ID))+
  geom_point(data = uMap_layout[is.na(term ) ], shape = 16, size=3, alpha=0.7, colour = 'grey70')+
  geom_point(data = uMap_layout[!is.na(term ) ],shape = 16, size=3, alpha=0.7, colour =colour_annot[8] )+
  
  geom_point(data = uMap_layout[Entry == to_check[4]],shape = 16, size=3, alpha=1, colour = colour_annot[4])+
  # geom_point(data = uMap_layout[!is.na(annot)][str_detect(annot,'^Int')], shape = 16, size=1, alpha=0.9)+
  # geom_point(data = uMap_layout[!is.na(annot)][str_detect(annot,'^Int',negate = T)], shape = 16, size=1, alpha=0.9)+
  theme(panel.background=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), legend.position = 'none',
        panel.border=element_rect(fill=NA, colour="black", size=0.95), panel.grid.major=element_blank(),
        axis.title=element_blank(), plot.margin = unit( c(0.5, 0.5, 0.5, 0.5), "mm"))+
  ggrepel::geom_text_repel(data = uMap_layout[!is.na(term ) ], colour  =colour_annot[8])+
  ggrepel::geom_label_repel(data = uMap_layout[Entry == to_check[4]], fill =colour_annot[4], colour = 'white' )+
  # scale_color_manual(values = c('grey50',colour_annot[4]))+
  # lims(y = c())
  scale_x_continuous(limits = c(umapx-0.3,umapx+0.3),expand = c(0,0))+
  scale_y_continuous(limits = c(umapy-0.3,umapy+0.3),expand = c(0,0))+
  ggtitle(glue::glue('{to_check[4]} {term_OI[2,Description]}'),
          subtitle = glue::glue('p.adj {format(term_OI[2,p.adjust],nsmall = 3)}')) 
ggsave(here::here('out','plots','umap_closeup_4.pdf'))


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
subset_stress = merge(subset_stress,Uniprot_mapping[,.(Entry,CGD)],by.x = 'Genes',by.y = 'CGD')
subset_stress= subset_stress[!is.na(Genes)]
subset_stress |> nrow() 
subset_stress[,Genes:=NULL]
subset_stress <- subset_stress |> tibble::column_to_rownames('Entry')
subset_stress = as.matrix(subset_stress)
tmp_medians <- apply( subset_stress, 1, median, na.rm = TRUE )  
subset_stress <- sweep( subset_stress, 1, tmp_medians, FUN = "-" )

Protein_OI = to_check[1]
position_term = 1
term_OI = enrichments_bicor[Protein == Protein_OI][Count>1][order(p.adjust)]
Interactors = term_OI[position_term,"geneID"] |> unlist() |> strsplit('/') |> unlist()
Interactors = Uniprot_mapping[ID %in% Interactors,Entry]
term = term_OI[position_term,"Description"] |> unlist()
go_Term_pval = term_OI[position_term,'p.adjust'] |> unlist()
bicorPO_present = subset_stress[rownames(subset_stress)==Protein_OI] |> is.finite()
bicorcor_OI = subset_stress[,bicorPO_present]
interactos_present = bicorcor_OI[str_detect(rownames(bicorcor_OI), paste(Interactors,collapse = '|')),] |> 
  is.finite() |> matrixStats::colSums2()
interactos_present = interactos_present>=quantile(interactos_present)[3]
bicorcor_OI_int = bicorcor_OI[,interactos_present]
bicorcor_OI_int  = bicorcor_OI_int[(bicorcor_OI_int |> 
                                          is.finite() |> matrixStats::rowSums2())>ncol(bicorcor_OI_int)/2,]
bicorcor_OI_int |> dim()
# tmp_medians <- apply( bicorcor_OI_int, 2, median, na.rm = TRUE )  
# bicorcor_OI_int <- sweep( bicorcor_OI_int, 2, tmp_medians, FUN = "-" )
# bicorcor_OI_int |> boxplot()
bicorcor_OI_int = bicorcor_OI_int |> as.data.table(keep.rownames = 'Uniprot')
bicorcor_OI_int_piv  = bicorcor_OI_int[,c(1,2:200)] |> 
  melt(id.vars = 'Uniprot', variable.name = 'experiment',value.name = 'LogRatio')
bicorcor_OI_int_piv[,category:=dplyr::case_when(
  str_detect(Uniprot,Protein_OI)~Protein_OI,
  str_detect(Uniprot,paste(Interactors,collapse = '|'))~term,
  TRUE~'other'
  
)]

colour_pallet = c('grey95','#ffcccb','darkred')
names(colour_pallet) = c('other',term,Protein_OI)
bicorcor_OI_int_piv |> ggplot(aes(x = experiment, y= LogRatio,colour= category,group=Uniprot))+
  # geom_line(data = bicorcor_OI_int_piv[category == 'other'],alpha = 0.05)+
  geom_line(data = bicorcor_OI_int_piv[category == term],alpha = 0.7)+
  geom_line(data = bicorcor_OI_int_piv[category == Protein_OI],alpha = 1)+theme_bw()+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_colour_manual(values = colour_pallet)+
  # annotate('text',label = , 
  #          x =100,  y = 4, angle = 0)+
  # # ylim(-5,5)+
  ggtitle(glue::glue('{Protein_OI} interactors {term} fdr/n adj.pval = {fifelse(round(go_Term_pval,3)==0,"<0.001",as.character(round(go_Term_pval,4)))}'))

