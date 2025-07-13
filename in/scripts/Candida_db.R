if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("AnnotationForge")
BiocManager::install("biomaRt")

library(AnnotationForge);library(clusterProfiler);library('biomaRt')
makeOrgPackageFromNCBI(version = "0.3",
                       author = "Savvas Kourtis <v1skourt@ed.ac.uk>",
                       maintainer = "Savvas Kourtis <v1skourt@ed.ac.uk>",
                       outputDir = ".",
                       tax_id = '237561',
                       genus = "Candida",
                       species = "albicans")


install.packages(here::here("org.Calbicans.eg.db"), repos=NULL,type="source")
DT = fread(here::here('out','datasets','whole_stressed_covariations.gz'))
library(org.Calbicans.eg.db)

# Uniprot mapping
Uniprot_annot = fread(here::here('in','datasets','uniprotkb_taxonomy_id_237561_2025_07_08.tsv.gz'))
Uniprot_annot[,CGD:=str_remove(CGD,';[:print:]*$')]
Uniprot_annot[,GeneID:=str_remove(GeneID,';[:print:]*$')]
covar_quantiles = DT[!is.na(corr_bicor),corr_bicor] |> quantile(probs = seq(0, 1, 0.001), type = 5)

covariation_partners <- DT[corr_bicor>=0.507,.(Protein_1,Protein_2,corr_bicor)]
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

gene2go = fread(here::here('gene2go.gz'))
gene2go[GeneID=='30514957']
crypto_214684 = gene2go[`#tax_id`==
          '214684']
unique_taxa = gene2go$`#tax_id` |> unique()
unique_taxa |> as.character() |> str_subset('214684')


enrichment_function <- function(Uniprot,Universe,covariations){
 i  = universe_geneID[Entry == Uniprot,GeneID]
 interactors = covariation_partners[GeneID.x == i|
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
    fwrite(term_OI, here::here('out','datasets',glue::glue('enrichment_terms_bicor.csv')), append = T)
    
  }
  
}
 }
}

plan(sequential)
plan(multisession, workers = 40)

furrr::future_walk(.x = universe_geneID$Entry, ~enrichment_function(Uniprot = .x,
                                                                    Universe = universe_geneID,
                                                                    covariations =covariation_partners ))


functions_candida = fread( here::here('out','datasets',glue::glue('enrichment_terms_bicor.csv')))
examples_to_show = merge(Uniprot_annot[,.(Entry,Annotation,`Protein existence`)],
      functions_candida[,.(min_padj = min(p.adjust),
                           N_terms = .N), by = Protein], by.x = 'Entry',
      by.y = 'Protein')
functions_candida$Protein |> unique()
