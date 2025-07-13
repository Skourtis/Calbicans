# identify c.candida Uniprot/gene IDs to build reactome/ string gold standards
library(data.table);library(stringr)
Sub_conditions = fread(here::here('in','datasets','subset_conditions_proteomes_refreplicates.tsv'))
# candida IDs
CGDs = unique(Sub_conditions$Genes)

# https://www.uniprot.org/uniprotkb?query=(taxonomy_id:237561)
Uniprot_annot = fread(here::here('in','datasets','uniprotkb_taxonomy_id_237561_2025_07_04.tsv.gz'))
Uniprot_annot[,CGD:=str_remove(CGD,';[:print:]*$')]
# 83 Identified Genes not in Uniprot
table(CGDs %in% as.vector(str_split(Uniprot_annot[CGD != '',CGD],';',simplify = T)))
Uniprots_detected = Uniprot_annot[CGD %in% CGD,Entry]

# STRING

STRING = fread(here::here('in','datasets','237561.protein.links.detailed.v12.0.txt.gz')) 
STRING[,`:=`(protein1 = str_remove(protein1,'237561.'),
             protein2 = str_remove(protein2,'237561.'))]
# It appears that all pairs are annotated as both A <-> B and B <-> A in this file. 
# Keep only one of these duplicates (could just to A > B but this is the safer option)
STRING[, Protein_1_sorted := ifelse(protein1 > protein2, protein1, protein2) ]           # Sort the IDs row-wise
STRING[, Protein_2_sorted := ifelse(protein1 < protein2, protein1, protein2) ]           # Sort the IDs row-wise
STRING <- STRING[, .(protein1 = Protein_1_sorted, protein2 = Protein_2_sorted, neighborhood, fusion, cooccurence,
                     coexpression, experimental, database, textmining, combined_score )]
STRING <- unique(STRING)     # Keep only unique pairs
# Restrict STRING table to pairs that have been in our dataset
STRING <- STRING[ protein1 %in% Uniprots_detected & protein2 %in% Uniprots_detected ]   


# Assign TP status to the existing pairs (will restrict to more significant interactions later)
STRING$STRING_Class <- "TP"

# To define FP pairs, find all combinations of these proteins that were not listed as associated
all_combinations <- data.table( t( combn( STRING[, unique(c(protein1, protein2))] , 2 )))
names(all_combinations) <- c("protein1", "protein2")
all_combinations[, Protein_1_sorted := ifelse(protein1 > protein2, protein1, protein2) ]            # Sort the IDs row-wise
all_combinations[, Protein_2_sorted := ifelse(protein1 < protein2, protein1, protein2) ]            # Sort the IDs row-wise
all_combinations <- all_combinations[, .(protein1 = Protein_1_sorted, protein2 = Protein_2_sorted)]

# Merge the tables, so that it now comprises all possible combinations of these proteins
STRING <- merge( STRING, all_combinations, by = c("protein1", "protein2"), all = TRUE )

# Assign FP pairs based on the fact that they are not TP pairs
STRING[ is.na(STRING_Class) , STRING_Class := "FP" ]

# The current definition of TP is probably still too loose, so I only include STRING interactions mapped with high confidence
# STRING interactions with low to medium confidence are completely excluded (i.e. they are also not used as false positives)
# Set low to medium confidence hits to NA
# Remove them from the table
# we decided to remove co-expre and leave threshold to 300

STRING_nocoexp = copy(STRING)
STRING_nocoexp = STRING_nocoexp[,combined_score_wo_coex:= combined_score-coexpression, by = .(protein1,protein2) 
]
STRING_nocoexp$combined_score_wo_coex |> hist()
STRING_nocoexp |> ggplot(aes(x = combined_score_wo_coex ))+
  geom_histogram( fill = 'grey90',colour = 'grey20')+
  geom_vline(colour = 'red', xintercept = 450)+
  theme_bw()+
  ggtitle('Defining a threshold for TPs without coexpression')
ggsave(here::here('out','plots','STRING_scores_GS definition.pdf'))
STRING_nocoexp = STRING_nocoexp[ combined_score < 450 , STRING_Class := NA 
][ !is.na( STRING_Class ) ] 
setnames(STRING_nocoexp,'STRING_Class','class')
fwrite(STRING_nocoexp[complete.cases(class),.(protein1,protein2,class)],here::here('out','datasets','STRING_TP_FP.gz') )

