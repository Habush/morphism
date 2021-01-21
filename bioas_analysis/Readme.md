Fuzzy-membership-based property vectors embedding (FMBPV)
---------------------------------------------------------

Patients data analysis 

- Requiremnts from OpenCog:

    Atomspace, PLN, URE, Agi-bio

    Note: Add the following line and rebuild Atomspace inorder for patterns like `(Not (set x))` to work

    https://github.com/tanksha/atomspace/commit/0ec8933d943842dccc04773f86273d5fcf299cfa

- Python dependencies

    scipy, sklearn, gensim, pandas, numpy

- Example:
Use the 10 patients sample dataset and generate an embedding vector

    `python patients_analysis.py --datapath cancer_data/sample_10_patients_data --outputpath /path/to/outputfile/` 

These will create the following in the outputpath directory given above
 - two csv files (before and after kpca) of property vector
 - AttractionLinks file
 - Zerovector file (only if there are patients with no attraction link to any property)
