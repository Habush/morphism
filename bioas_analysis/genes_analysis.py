import os
from utils import *
from generate_embedding import *
from preprocess import calculate_truth_values, infer_subsets
import argparse

log.set_level("ERROR")

def load_atomspace(datapath, atomspace):
  print("--- Load preprocessed data to the AtomSpace")
  for i in os.listdir(datapath):
    if i.endswith(".scm"):
      scheme_eval(atomspace, "(load-file \"{}\")".format(os.path.join(datapath, i)))  

def generate_subsets(atomspace):
  print("--- Infer subsets between gene singletons and concepts")
  scheme_eval(atomspace, "(load-from-path \"rules/subset-genes-rule.scm\")")
  scheme_eval(atomspace, "(create-subset-lns 'GeneNode)")

def gene_singleton_tv(atomspace):
    print("--- Assign TV for gene singletons")
    all_genes = atomspace.get_atoms_by_type(types.GeneNode)
    universe_size = len(all_genes)
    tv_confidence = float(scheme_eval(atomspace, "(count->confidence " + str(universe_size) + ")"))
    strength = 1/universe_size
    scheme_eval(atomspace, 
    "(map (lambda (a) (cog-set-tv! (Set a) (stv {} {}))) (cog-get-atoms 'Gene))".format(strength, tv_confidence))

def infer_negation(atomspace):
  print("--- Inferring Subset Negation")
  # (Subset (Set A) B) |- (Subset (Not (Set A) B))
  scheme_eval(atomspace, "(load-from-path \"rules/subset_negation_rule.scm\")")
  scheme_eval(atomspace, "(create-subset-neg-lns 'GeneNode)")

def infer_attractions(atomspace):
  print("--- Inferring AttractionLinks")
  # (Subset (Set A) B)) (Subset (Not (Set A) B)) |- (Attraction (Set A) B)
    #:complexity-penalty 10)"""
  scheme_eval(atomspace, "(load-from-path \"rules/subset_attraction_rule.scm\")")
  scheme_eval(atomspace, "(create-attr-lns 'GeneNode)")

def get_intentional_similarity(atomspace):
  # (Attraction (Set A) X)) (Attraction (Set B) X)) |- Intensional-similarity (Set A) (Set B))
  scheme_eval(atomspace, "(pln-load 'empty)")
  scheme_eval(atomspace, "(load-from-path \"rules/genes_intensional_similarity.scm\")")
  scheme_eval(atomspace, "(create-ints-similarity-lns)")

def export_result(datapath, atomspace):
  print("--- {} Exporting Atoms to files")
  result_scm = os.path.join(datapath, "AttractionLinks-results.scm.bc")
  att = atomspace.get_atoms_by_type(types.AttractionLink) + atomspace.get_atoms_by_type(types.IntensionalSimilarityLink)
  with open(result_scm, "w") as f:
    f.write("\n".join([str(a) for a in att]))

def generate_links_genes(datapath, doinhrule=False, genes=False):
  ### Initialize the AtomSpace ###
  atomspace = AtomSpace()
  initialize_opencog(atomspace)

  scheme_eval(atomspace, "(add-to-load-path \".\")")
  scheme_eval(atomspace, "(use-modules (opencog) (opencog exec) (opencog logger) (opencog bioscience) (opencog ure) (opencog pln) (opencog persist-file) (srfi srfi-1))")
  # scheme_eval(atomspace, """(ure-logger-set-level! "debug")""")
  load_atomspace(datapath, atomspace)
  if genes:
    for g in atomspace.get_atoms_by_type(types.GeneNode):
      if not g.name in genes:
        scheme_eval(atomspace,"(cog-extract-recursive! {})".format(g))
  if doinhrule:
    infer_subsets(atomspace)
  calculate_truth_values(atomspace) # for the memberlink and GO
  gene_singleton_tv(atomspace)
  generate_subsets(atomspace)
  infer_negation(atomspace)
  infer_attractions(atomspace)
  get_intentional_similarity(atomspace)
  return atomspace

def parse_args():
    parser = argparse.ArgumentParser(description='Generate embedding vector for patients')
    parser.add_argument('--datapath', type=str, default='',
                        help='path to the background knowlegde data for gene (e.g membership to GO/Pathway)')
    parser.add_argument('--outputpath', type=str, default='',
                        help='path to store the output files')
    parser.add_argument('--doinhrule', type=bool, default=False,
                        help='Apply inheritance rule and Propagate the memberlinks relationship to the parent GO')
    parser.add_argument('--genes', type=str, default='',
                        help='Text file with list of gene names to generate embedding for')
    return parser.parse_args()

if __name__=="__main__":
  arguments = parse_args()
  if arguments.datapath:
    base_datasets_dir = arguments.datapath
  else:
    base_datasets_dir = os.getcwd() + "/cancer_data/"
  if arguments.outputpath:
    base_results_dir = arguments.outputpath
  else:
    base_results_dir = os.getcwd() + "/results/cancer-{}/".format(str(date.today()))
  if arguments.doinhrule:
    doinhrule = True
  else:
    doinhrule = False
  if arguments.genes:
    genes = open(arguments.genes,"r").read().splitlines()
  else:
    genes = False
  kbs = generate_links_genes(base_datasets_dir, doinhrule=doinhrule, genes=genes)
  export_result(base_results_dir, kbs)
  generate_embeddings(base_results_dir,"GeneNode", kb_atomspace=kbs)
