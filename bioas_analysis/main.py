from preprocess import *
from concepts_analysis import *
from genes_analysis import *
from generate_embedding import *

output_path = generate_atoms()
generate_attractionlinks_concepts(output_path)
generate_attractionlinks_genes(output_path)
generate_embeddings("FMBPV",output_path)
