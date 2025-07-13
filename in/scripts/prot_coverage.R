# protein-coverage proteomehd2
library(data.table);library(stringr);library(ggplot2)

# loading fasta file for sequences
isoform_fasta <- seqinr::read.fasta('/home/v1skourt/proteomehd2/fasta/micro_prot_decoy.fasta',
                                    seqtype = "AA",as.string = TRUE)
isoform_fasta = data.table(Annotation  = unlist(seqinr::getAnnot(isoform_fasta)) |>  str_remove(' [:print:]*$') |> str_remove('>'),
                           Seq = seqinr::getSequence(isoform_fasta, as.string = T) |>  unlist())
fasta_pgFDR = copy(isoform_fasta)
fasta_pgFDR[,pro_length:= str_length(Seq)]

# loading pgFDR for proteins which survive FDR
pgFDR_pGroups =fread(here::here('in','datasets',
                              'picked_protein_group_no_remap',
                              'proteinGroups_picked_protein_group_no_remap.txt'))
pgFDR_pGroups[,Annotation:= str_remove(`Majority protein IDs`,';[:print:]*$')]
fasta_pgFDR = merge(fasta_pgFDR,pgFDR_pGroups, by = 'Annotation', all.x = F)
fasta_pgFDR[,surviving_FDR:= factor(fifelse(`Q-value`<0.01,'pgFDR<0.01','pgFDR>0.01'),
                                    levels = c('pgFDR>0.01','pgFDR<0.01'))]
breaks_tmp = seq(0,400,by = 50)

# defining protein size bins to check pgFDR across protein size
fasta_pgFDR[,binned_size:= cut(pro_length,breaks = breaks_tmp)]
fasta_pgFDR[,amino_acid_length:= as.character(max(pro_length)), by =binned_size  ]
fasta_pgFDR[,type:= fifelse(str_detect(Annotation,'^(REV_sp|sp)'),'swissprot','microprotein')]

# average per bin
perc_surv = fasta_pgFDR[,.(total_per_class = .N),by = .(type,amino_acid_length,surviving_FDR)]
perc_surv[,total:= sum(total_per_class ), by = .(type,amino_acid_length)]
perc_surv[,perc_surv :=total_per_class/total ]
order_bins = fasta_pgFDR[order(pro_length)][,head(.SD,1), by = binned_size][,amino_acid_length]

surviving_50aa_micro = fasta_pgFDR[!is.na(`Q-value`)][Reverse != '+' ][amino_acid_length %in% c('100')] |> 
  ggplot(aes(x = type , fill =surviving_FDR ))+
  geom_bar()+
  # facet_wrap('type',nrow = 2, scales = 'free_y')+
  geom_text(data = perc_surv[surviving_FDR  =='pgFDR<0.01'][amino_acid_length %in% c('100')],size= 1.8,
            y = -25,aes( label =paste0(round(perc_surv*100,0),'%' )))+
  scale_fill_manual(values =c('grey80','grey30'))+
  labs(x = 'Prots with length between 50-100 aa', y  = 'Number of proteins')

ProHD2_figs[['surviving_50aa_micro']] <- surviving_50aa_micro
ProHD2_figdata[['surviving_50aa_micro']] <- fasta_pgFDR[!is.na(`Q-value`)
                                                     ][Reverse != '+' ][amino_acid_length %in% c('100')]
all_pgFDR_prots_size = fasta_pgFDR[!is.na(`Q-value`)][Reverse != '+' ] |> 
  ggplot(aes(x = factor(amino_acid_length,levels =order_bins)  ))+
  geom_bar()+
  facet_wrap('type',nrow = 2, scales = 'free_y')+
  # geom_text(data = perc_surv[surviving_FDR  ==T],y = -100,aes( label =paste0(round(perc_surv*100,0),'%' )))+
  labs(x = 'Max prot length per bin', y  = 'Number of proteins')
ProHD2_figs[['all_pgFDR_prots_size']] <- all_pgFDR_prots_size
ProHD2_figdata[['all_pgFDR_prots_size']] <- fasta_pgFDR[!is.na(`Q-value`)][Reverse != '+' ]

pgFDR_swissprot_bias = fasta_pgFDR[!is.na(`Q-value`)][Reverse != '+' ][type =='swissprot'] |> 
  ggplot(aes(x = factor(amino_acid_length,levels =order_bins) , fill =surviving_FDR ))+
  geom_bar(position = 'fill')+
  # facet_wrap('type',nrow = 2, scales = 'free_y')+
  geom_text(data = perc_surv[surviving_FDR  =='pgFDR<0.01'][type =='swissprot'],size= 1.8,
            y = -0.01,aes( label =paste0(round(perc_surv*100,0),'%' )))+
  labs(x = 'Max prot length per bin', y  = 'Fraction of swissprot proteins')+
  scale_fill_manual(values =c('grey80','grey30'))
ProHD2_figs[['pgFDR_swissprot_bias']] <- pgFDR_swissprot_bias
ProHD2_figdata[['pgFDR_swissprot_bias']] <-fasta_pgFDR[!is.na(`Q-value`)][Reverse != '+' ][type =='swissprot']
ProHD2_figdata$order_bins = order_bins
# calculating protein coverage per protein
PG_all_cont = pgFDR_pGroups$`Protein IDs` |> strsplit(';') |> unlist() |> stringr::str_subset('orf|ribo|ORF|cogn',negate = T) |> stringr::str_remove('[:print:]*_') |> 
  table() |> names()
PG_all_cont = paste0('_',PG_all_cont[stringr::str_detect(PG_all_cont,'HUMAN',negate = T)])
pgFDR_pGroups = pgFDR_pGroups[,`Potential contaminant`:= str_detect(`Protein IDs`,paste(PG_all_cont,collapse = '|'))]
pgFDR_pGroups = pgFDR_pGroups[`Q-value`<0.01]
pgFDR_pGroups = pgFDR_pGroups[Reverse != '+' ]
pgFDR_pGroups = pgFDR_pGroups[`Potential contaminant`==F]

pgFDR_pGroups = str_split(pgFDR_pGroups[,`Protein IDs`] ,';')

length_PG_per_group = purrr::map(.x = pgFDR_pGroups,~
                                   purrr::map_dbl(.x= .x,~isoform_fasta[Annotation==.x,Seq] |> str_length()))
longest_PG = purrr::map_dbl(.x = length_PG_per_group,~head(which(.x ==max(.x)),1))
longest_PGs = purrr::map2_chr(.x = pgFDR_pGroups,.y = longest_PG,~.x[.y])

Protein_coverage = isoform_fasta[Annotation %in% longest_PGs]

psm_percolator = fread(here::here('in','datasets','pep_to_prot_mapping.txt'),header = F)
names(psm_percolator) = c('sequence','Uniprot')
psm_percolator[,sequence:= str_remove(sequence,'^n')]
setkey(psm_percolator,'Uniprot')

count_coverage <- function(Annotation,sequence,psms_percolator){
  
  # Annotation = Protein_coverage[2,1] |> unlist()
  Annotation = Annotation |> str_remove('^[:print:]*?\\|') |> str_remove('\\|[:print:]*') 
  # sequence  = Protein_coverage[2,2] |> unlist()
  psm_percolator_tmp = psm_percolator[str_detect(Uniprot,glue::glue('(^|;){Annotation}')),sequence]
  positions_tmp = c()
  # print(Annotation)
  # print(sequence)
  for(i in psm_percolator_tmp){
    # sometimes it's missing for leucine and isoleucine
    # i = psm_percolator_tmp[1]
    location = str_locate(sequence,i)
    if(!any(is.na(location))){
      positions_tmp = c(positions_tmp,location[1]:location[2])}
  }
  positions_tmp = unique(positions_tmp)
  perc_coverage = length(positions_tmp)/str_length(sequence)
  perc_coverage
}
Protein_coverage[,perc_coverage:= count_coverage(Annotation,Seq,psms_percolator),by = Annotation]
fwrite(Protein_coverage,here::here('out','datasets','pgFDR_prot_coverage_.tsv' ))

ggplot(Protein_coverage,aes(x = perc_coverage))+
  geom_histogram(bins = 50)+
  ggtitle('protein coverage prohd2')

Protein_coverage = fread(here::here('out','datasets','pgFDR_prot_coverage_.tsv' ))
Protein_coverage_plot = Protein_coverage|> ggplot(aes(x  =perc_coverage))+
  geom_histogram(bins = 50, fill = 'grey50')+
  geom_vline(xintercept = median(Protein_coverage$perc_coverage,na.rm = T), linetype = 'dashed')+
  labs(x = 'Proteins sequence coverage', y ='Number of Proteins')
ggsave(here::here('out','plots','protein_coverage.pdf'))

ProHD2_figdata[['Protein_coverage']] <- Protein_coverage
ProHD2_figs[['Protein_coverage']] <- Protein_coverage_plot


