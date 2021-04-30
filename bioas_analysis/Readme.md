Fuzzy-membership-based property vectors embedding (FMBPV)
---------------------------------------------------------

- Requiremnts from OpenCog:

    Atomspace, PLN, URE, Agi-bio

    Note: Add the following line and rebuild Atomspace inorder for patterns like `(Not (set x))` to work

    https://github.com/tanksha/atomspace/commit/0ec8933d943842dccc04773f86273d5fcf299cfa

- Python dependencies

    scipy, sklearn, gensim, pandas, numpy

Patients data analysis 
- The python script dopln_generate_embedding.py does the following
    - accepts the RAW gene expression dataset (datapath)
    - Normalize the dataset (to scale the expression values to between 0 and 1 tobe used as STV)
    - Converts it into Atomese 
    - Apply different PLN rules (if you want to infer GO and pathways set the infergo parameter to True)
    - Generate embedding vectors for Patients (if you want to generate embedding vectors from inferred GO and pathways only, set plnonly to True)

Genes Analysis
- The python script genes_analysis.py does the following
    - Read the scm files (background knowledge for Gene to GO/pathway memberlinks)
    - Apply PLN rule to infer AttractionLinks between Genes and GO/pathway, and Intensional similarity between genes
      Note, if you want to apply inheritance rule and propagate the Membership relationships to Parent GO and pathways, set doinhrule to True
    - Generate embedding vectors for Genes

check each files for other optional parameters to set

