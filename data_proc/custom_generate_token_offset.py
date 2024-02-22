# 生成offset文件和tokens文件
import torch
import pickle
token_dim = 5120
species_to_offsets = {}

all_pe = torch.load("/kaggle/working/all_tokens.torch")[0:4] # read in existing token file to make sure 
# that special vocab tokens are the same for different seeds

offset = len(all_pe) # special tokens at the top!

PE = torch.load(SPECIES_PROTEIN_EMBEDDINGS_PATH)

pe_stacked = torch.stack(list(PE.values()))
print(all_pe.shape)
print(pe_stacked.shape)
all_pe = torch.vstack((all_pe, pe_stacked))

species_to_offsets[species] = offset

print("CHROM_TOKEN_OFFSET:", all_pe.shape[0])
torch.manual_seed(TAXONOMY_ID)
CHROM_TENSORS = torch.normal(mean=0, std=1, size=(N_UNIQ_CHROM, 5120)) 
# N_UNIQ_CHROM is the total number of chromosome choices, it is hardcoded for now (for species in the training data)
all_pe = torch.vstack(
    (all_pe, CHROM_TENSORS))  # Add the chrom tensors to the end
all_pe.requires_grad = False


torch.save(all_pe, f"/kaggle/working/{SPECIES_NAME}_pe_tokens.torch")

with open(f"/kaggle/working/{SPECIES_NAME}_offsets.pkl", "wb+") as f:
    pickle.dump(species_to_offsets, f)
print("Saved PE, offsets file")
print(f"/kaggle/working/{SPECIES_NAME}_offsets.pkl")
print(SPECIES_PROTEIN_EMBEDDINGS_PATH)
