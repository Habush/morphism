import os
from datetime import date, datetime
from utils import *
from generate_embedding import *
import argparse

def populate_atomspace(atomspace, path):
  print("--- {} Populating the AtomSpace {}".format(datetime.now(), path))
  for i in os.listdir(path):
    if str(i).endswith(".scm"):
      print(i)
      scheme_eval(atomspace, "(load-file \"{}\")".format(os.path.join(path, i)))  

def preprocess(atomspace):
  print("--- {} Preprocessing the AtomSpace".format(datetime.now()))
  # remove gene expressions, pln results and patient status with mean 0
  for e in atomspace.get_atoms_by_type(types.SubsetLink) + atomspace.get_atoms_by_type(types.EvaluationLink):
    if e.tv.mean == 0:
      scheme_eval(atomspace,"(cog-delete {})".format(e))

  # remove patients with no gene expression
  filter_out = []
  all_patients = atomspace.get_atoms_by_type(types.PatientNode)
  with_geneexpression = set([x.out[1] for x in atomspace.get_atoms_by_type(types.EvaluationLink) if x.out[0].type == types.LazyExecutionOutputLink])
  for p in all_patients:
    if p not in with_geneexpression:
      filter_out.append(str(p))
      scheme_eval(atomspace, "(cog-delete-recursive! {})".format(p))
  with open(base_results_dir + "filtered_out_patients.scm", "w") as f:
    f.write("\n".join(filter_out))
  
  return len(with_geneexpression)

def apply_subset_rule1(atomspace):
  print("--- {} Inferring subsets, rule 1".format(datetime.now()))
  scheme_eval(atomspace, "(pln-load 'empty)")
  scheme_eval(atomspace, "(pln-load-from-path \"rules/patients_subset_rule.scm\")")
  scheme_eval(atomspace, "(pln-add-rule \"patient-data-subset-rule\")")
  target = """
    (Subset 
      (Set (Variable "$p"))
      (Variable "$ppty"))"""

  bc = """
    (pln-bc 
      {}
      #:vardecl (VariableSet 
        (TypedVariable (Variable "$p") (Type "PatientNode")) 
        (TypedVariable (Variable "$ppty") (Type "SatisfyingSetScopeLink")))
      #:maximum-iterations 2
      #:complexity-penalty 10)
      """.format(target)

  scheme_eval(atomspace, bc) 
  scheme_eval(atomspace, "(pln-load 'empty)")
  scheme_eval(atomspace, "(pln-load-from-path \"rules/patients_subset_rule.scm\")")
  scheme_eval(atomspace, "(pln-add-rule \"patient-data-boolean-subset-rule\")")
  scheme_eval(atomspace, bc) 

  scheme_eval(atomspace, "(pln-load 'empty)")
  scheme_eval(atomspace, "(pln-load-from-path \"rules/patients_subset_rule.scm\")")
  scheme_eval(atomspace, "(pln-add-rule \"gene-expression-subset-rule\")")
  scheme_eval(atomspace, bc) 
  
def apply_subset_rule2(atomspace):
  print("--- {} Inferring subsets, rule 2".format(datetime.now()))
  scheme_eval(atomspace, "(pln-load 'empty)")
  scheme_eval(atomspace, "(pln-load-from-path \"rules/patients-ppty-rule.scm\")")
  scheme_eval(atomspace, "(pln-add-rule \"patient-ppty-subset-rule\")")
  target = """
    (Subset 
      (Set (Variable "$p"))
      (Variable "$ppty"))"""

  bc = """
  (pln-bc 
    {}
    #:vardecl (VariableSet 
      (TypedVariable (Variable "$p") (Type "PatientNode")) 
      (TypedVariable (Variable "$ppty") (Type "SatisfyingSetScopeLink")))
    #:maximum-iterations 2
    #:complexity-penalty 10)""".format(target)
  scheme_eval(atomspace, bc)  

def generate_attraction_links(atomspace):
  print("--- {} Generate Attractions".format(datetime.now()))
  scheme_eval(atomspace, "(pln-load 'empty)")
  scheme_eval(atomspace, "(pln-load-from-path \"rules/subset_negation_rule.scm\")")
  scheme_eval(atomspace, "(pln-load-from-path \"rules/subset_attraction_rule.scm\")")
  scheme_eval(atomspace, "(pln-add-rule \"subset-negation-patients-rule\")")
  scheme_eval(atomspace, "(pln-add-rule \"subset-attraction-patients-rule\")")
  target = """
  (Subset
  (Not (Set (Variable "$X")))
  (Variable "$Y"))
  """
  target2 = """ 
  (Attraction 
      (Set (Variable "$X")) 
      (Variable "$Y"))
  """
  bc = """
  (pln-bc 
    {}

    #:maximum-iterations 2
    #:complexity-penalty 10)"""
  scheme_eval(atomspace, bc.format(target))
  scheme_eval(atomspace, bc.format(target2))

def remove_processed_subsets(atomspace):
  for e in atomspace.get_atoms_by_type(types.SubsetLink):
    if e.out[0].type != types.SetLink:
      scheme_eval(atomspace,"(cog-delete {})".format(e))

def export_all_atoms(atomspace, base_results_dir):
  print("--- {} Exporting Atoms to files".format(datetime.now()))
  result_scm = base_results_dir + "AttractionLinks-results.scm"
  att = atomspace.get_atoms_by_type(types.AttractionLink)
  with open(result_scm, "w") as f:
    f.write("\n".join([str(a) for a in att]))

def generate_atoms(base_results_dir, base_datasets_dir, filterbp):
    ### Initialize the AtomSpace ###
    atomspace = AtomSpace()
    initialize_opencog(atomspace)

    ### Guile setup ###
    scheme_eval(atomspace, "(add-to-load-path \".\")")
    scheme_eval(atomspace, """
    (use-modules (opencog) (opencog bioscience) (opencog ure)
    (opencog pln) (opencog persist-file) (srfi srfi-1) (opencog exec))
    """)

    populate_atomspace(atomspace,base_datasets_dir)
    total_patients = preprocess(atomspace)
    if filterbp:
      print("--- {} Filtering BP".format(datetime.now()))
      atomspace = filter_bp(atomspace)
    scheme_eval(atomspace, "(define total_patients {})".format(total_patients))
    apply_subset_rule1(atomspace)
    apply_subset_rule2(atomspace)
    remove_processed_subsets(atomspace)
    generate_attraction_links(atomspace)
    return atomspace

def parse_args():
    parser = argparse.ArgumentParser(description='Generate embedding vector for patients')
    parser.add_argument('--datapath', type=str, default='',
                        help='path to the source data')
    parser.add_argument('--outputpath', type=str, default='',
                        help='path to store the output files')
    parser.add_argument('--filterbp', type=str, default=False,
                        help='filter biological process only')
    return parser.parse_args()

if __name__ == "__main__":
  arguments = parse_args()
  if arguments.datapath:
    base_datasets_dir = arguments.datapath
  else:
    base_datasets_dir = os.getcwd() + "/cancer_data/"
  if arguments.outputpath:
    base_results_dir = arguments.outputpath
  else:
    base_results_dir = os.getcwd() + "/results/cancer-{}/".format(str(date.today()))
  if arguments.filterbp:
    filterbp = True
  else:
    filterbp = False

  if not os.path.exists(base_results_dir):
    os.makedirs(base_results_dir) 
  kb_as = generate_atoms(base_results_dir, base_datasets_dir, filterbp)
  generate_embeddings("FMBPV",base_results_dir,"PatientNode", kb_atomspace=kb_as)
  export_all_atoms(kb_as, base_results_dir)
  print("Done {}".format(datetime.now()))