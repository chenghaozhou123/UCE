# 进行测试时部分数据的生成 开始执行
import numpy as np
import pickle as pkl
import pandas as pd
SPECIES_NAME = "arabidopsis" # short hand name for this species, will be used in arguments and files

# Path to the species proteome
SPECIES_PROTEIN_FASTA_PATH = "/kaggle/working/Arabidopsis_thaliana.TAIR10.pep.all.fa"

# Path to the ESM2 Embeddings
SPECIES_PROTEIN_EMBEDDINGS_PATH = "/kaggle/input/pretrain-embedding/Arabidopsis_thaliana.TAIR10.pep.all.gene_symbol_to_embedding_ESM2.pt"

# primary_assembly name, this needs to be matched to the FASTA file
ASSEMBLY_NAME = "TAIR10"
# NCBI Taxonomy ID, please set this so that if someone else also embeds the same species,
# randomly generated chromosome tokens will be the same
TAXONOMY_ID = 3702  # 在官网的id
species_to_paths = {
    SPECIES_NAME: SPECIES_PROTEIN_FASTA_PATH,
}

species_to_ids = {
    SPECIES_NAME: ASSEMBLY_NAME,
}
all_pos_def = []

missing_genes = {}
for species in species_to_ids.keys():
    missing_genes[species] = []
    proteome_path = species_to_paths[species]
    species_id = species_to_ids[species]

    with open(proteome_path) as f:
        proteome_lines = f.readlines()

    gene_symbol_to_location = {}
    gene_symbol_to_chrom = {}

    for line in proteome_lines:
        if line.startswith(">"):
            split_line = line.split()
            gene_symbol = [token for token in split_line if token.startswith("gene_symbol")]
            if len(gene_symbol) > 0:
                gene_symbol = gene_symbol[0].split(":")
                
                if len(gene_symbol) == 2:
                    gene_symbol = gene_symbol[1]
                elif len(gene_symbol) > 2:
                    gene_symbol = ":".join(gene_symbol[1:]) # fix for annoying zebrafish gene names with colons in them
                else:
                    1/0 # something weird happening, throw an error
                
                
                chrom = None
                
                chrom_arr = [token for token in split_line if token.startswith("chromosome:")]
                if len(chrom_arr) > 0:
                    chrom = chrom_arr[0].replace("chromosome:", "")
                else:
                    chrom_arr = [token for token in split_line if token.startswith("primary_assembly:")]
                    if len(chrom_arr) > 0:
                        chrom = chrom_arr[0].replace("primary_assembly:", "")
                    else:
                        chrom_arr = [token for token in split_line if token.startswith("scaffold:")] 
                        if len(chrom_arr) > 0:
                            chrom = chrom_arr[0].replace("scaffold:", "")
                if chrom is not None:
                    gene_symbol_to_location[gene_symbol] = chrom.split(":")[2]
                    gene_symbol_to_chrom[gene_symbol] = chrom.split(":")[1]
                else:
                    missing_genes[species].append(gene_symbol)
                    

    positional_df = pd.DataFrame()
    positional_df["gene_symbol"] = [gn.upper() for gn in list(gene_symbol_to_chrom.keys())]
    positional_df["chromosome"] = list(gene_symbol_to_chrom.values())
    positional_df["start"] = list(gene_symbol_to_location.values())
    positional_df = positional_df.sort_values(["chromosome", "start"])
    #positional_df = positional_df.set_index("gene_symbol")
    positional_df["species"] = species
    all_pos_def.append(positional_df)
    
master_pos_def = pd.concat(all_pos_def)
print(master_pos_def)
master_pos_def["species"].value_counts() # double check how many genes are mapped



for k, v in missing_genes.items():
    print(f"{k}: {len(v)}") # are any genes missing?


# Count genes per chromosome 保存文件
for species in species_to_ids.keys():
    print("*********")
    print(species)
    display(master_pos_def[master_pos_def["species"] == species]["chromosome"].value_counts().head(50))
    print("*********")
    
master_pos_def.to_csv(f"/kaggle/working/{SPECIES_NAME}_to_chrom_pos.csv", index=False) # Save the DF
# The chromosome file path will be:
print(f"{SPECIES_NAME}_to_chrom_pos.csv")
N_UNIQ_CHROM = len(master_pos_def[master_pos_def["species"] == species]["chromosome"].unique())
print(N_UNIQ_CHROM)
