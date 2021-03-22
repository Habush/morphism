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

def preprocess(atomspace, universe):
  print("--- {} Preprocessing the AtomSpace".format(datetime.now()))
  # remove gene expressions, pln results and patient status with mean 0
  for e in atomspace.get_atoms_by_type(types.SubsetLink) + atomspace.get_atoms_by_type(types.EvaluationLink):
    if e.tv.mean == 0:
      scheme_eval(atomspace,"(cog-delete {})".format(e))

  all_patients = atomspace.get_atoms_by_type(types.PatientNode)  
  return len(all_patients)

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
  result_scm = base_results_dir + "AttractionLinks-results.scm.bc"
  att = atomspace.get_atoms_by_type(types.AttractionLink)
  with open(result_scm, "w") as f:
    f.write("\n".join([str(a) for a in att]))

def generate_atoms(base_results_dir, base_datasets_dir, filterbp, universe=False):
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
    total_patients = preprocess(atomspace, universe)
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
    parser.add_argument('--universe', type=str, default='',
                        help='Text file containing list of patients to do embedding for (train or test sets)')
    parser.add_argument('--norm', type=str, default='',
                        help='Type of normalization none, scale, standard or quantile')
    parser.add_argument('--p', type=str, default='',
                        help='The value for p in the weight product s^p * c^(T-p)')
    parser.add_argument('--T', type=str, default='',
                        help='The value for T in the weight product s^p * c^(T-p)')
    parser.add_argument('--vec_space', type=str, default='',
                        help='The initial vector space to use as a base')
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
  if arguments.universe:
    universe_patients = open(arguments.universe, "r").read().splitlines()
  else:
    universe_patients = False 
  if arguments.norm:
    norm = arguments.norm
  else:
    norm = False 
  if arguments.p:
    p = float(arguments.p)
  else:
    p = 1
  if arguments.T:
    T = float(arguments.T)
  else:
    T = 2
  if arguments.vec_space:
    vec_space = pd.read_csv(arguments.vec_space, sep="\t")
    test_dataset = True
  else:
    test_dataset = False

  if not os.path.exists(base_results_dir):
    os.makedirs(base_results_dir)

  kb_as = generate_atoms(base_results_dir, base_datasets_dir, filterbp, universe=universe_patients)
  export_all_atoms(kb_as, base_results_dir) 

  if test_dataset:
    generate_embeddings(base_results_dir,"PatientNode", kb_atomspace=kb_as, p=p, T=T, normalization=norm, 
                        test_data=True, vector_space=vec_space)  
  else:
    generate_embeddings(base_results_dir,"PatientNode", kb_atomspace=kb_as, p=p, T=T, normalization=norm)
  print("Done {}".format(datetime.now()))