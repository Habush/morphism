import os
from datetime import date
from utils import *
from generate_embedding import *
import argparse

def populate_atomspace(atomspace, path):
  print("--- Populating the AtomSpace {}".format(path))
  for i in os.listdir(path):
    if str(i).endswith(".scm"):
      print(i)
      scheme_eval(atomspace, "(load-file \"{}\")".format(os.path.join(path, i)))  

def preprocess(atomspace):
  # remove gene expressions and patient status with stv 0
  for e in atomspace.get_atoms_by_type(types.EvaluationLink):
    if e.tv.mean == 0:
      scheme_eval(atomspace,"(cog-delete {})".format(e))

  # remove patients with no gene expression
  filter_out = []
  all_patients = [x for x in atomspace.get_atoms_by_type(types.ConceptNode) if str(x.name).isnumeric() and len(str(x.name)) > 1]
  with_geneexpression = set([x.out[1] for x in atomspace.get_atoms_by_type(types.EvaluationLink) if x.out[0].type == types.LazyExecutionOutputLink])
  for p in all_patients:
    if p not in with_geneexpression:
      filter_out.append(str(p))
      scheme_eval(atomspace, "(cog-delete-recursive! {})".format(p))
  with open(base_results_dir + "filtered_out_patients.scm", "w") as f:
    f.write("\n".join(filter_out))
  
  return len(with_geneexpression)

def apply_subset_rule1(atomspace):
  print("--- Inferring subsets")
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
        (TypedVariable (Variable "$p") (Type "ConceptNode")) 
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
  scheme_eval(atomspace, "(pln-load 'empty)")
  scheme_eval(atomspace, "(pln-load-from-path \"rules/patients-ppty-rule.scm\")")
  scheme_eval(atomspace, "(pln-add-rule \"patient-ppty-subset-rule\")")
  target = """ 
  (Subset 
      (Set (Variable "$p"))
      (SatisfyingSet
          (Variable "$pt")
          (SubsetLink
          (And
              (Variable "$c")
              (ConceptNode "profiled-genes"))
          (SatisfyingSetScopeLink
              (VariableNode "$G")
              (EvaluationLink
              (LazyExecutionOutputLink
                  (Variable "$S")
                  (VariableNode "$G"))
              (Variable "$pt"))))))"""

  bc = """
  (pln-bc 
    {}
    #:maximum-iterations 2
    #:complexity-penalty 10)""".format(target)
  scheme_eval(atomspace, bc)  

def generate_attraction_links(atomspace):
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

def calculate_truth_values(atomspace, len_patients):
  print("--- Calculating Truth Values")

  def get_confidence(count):
    return float(scheme_eval(atomspace, "(count->confidence {})".format(str(count))))

  total_patients = len_patients
  ppty_subset = [s for s in atomspace.get_atoms_by_type(types.SubsetLink) 
      if s.out[0].type == types.SetLink and s.out[1].type == types.SatisfyingSetScopeLink]
  for s in ppty_subset:
    try:
        # set tv of patients as 1 / total number of patients, 
        strength = 1 / total_patients
        confidence = get_confidence(total_patients)
        s.out[0].tv = TruthValue(strength, confidence)

        # get the SatisfyingSetScopeLink and use Get to get the list, devide the number of elements by total no of patients 
        satisf = s.out[1]
        outgoing = " ".join([str(i) for i in satisf.out])
        strength = int(scheme_eval(atomspace, """(length (cog-outgoing-set (cog-execute! (Get {}))))""".format(outgoing))) / total_patients
        confidence = get_confidence(total_patients)
        s.out[1].tv = TruthValue(strength, confidence)
    except Exception as e:
        print(e)
        continue

def remove_processed_subsets(atomspace):
  for e in atomspace.get_atoms_by_type(types.SubsetLink):
    if e.out[0].type != types.SetLink:
      scheme_eval(atomspace,"(cog-delete {})".format(e))

def export_all_atoms(atomspace, base_results_dir):
  print("--- Exporting Atoms to files")
  subset_links_scm = base_results_dir + "subset-links.scm"
  attraction_links_scm = base_results_dir + "attraction-links.scm"
  write_atoms_to_file(subset_links_scm, "(cog-get-atoms 'SubsetLink)", atomspace)
  write_atoms_to_file(attraction_links_scm, "(cog-get-atoms 'AttractionLink)", atomspace)

def generate_atoms(base_results_dir, base_datasets_dir):
    ### Initialize the AtomSpace ###
    atomspace = AtomSpace()
    initialize_opencog(atomspace)

    ### Guile setup ###
    scheme_eval(atomspace, "(add-to-load-path \".\")")
    scheme_eval(atomspace, """
    (use-modules (opencog) (opencog bioscience) (opencog ure) (opencog logger)
    (opencog pln) (opencog persist-file) (srfi srfi-1) (opencog exec))
    (ure-logger-set-level! "debug")
    """)
    scheme_eval(atomspace, " ".join([
    "(define (write-atoms-to-file file atoms)",
        "(define fp (open-output-file file))",
        "(for-each",
        "(lambda (x) (display x fp))",
        "atoms)",
        "(close-port fp))"]))

    populate_atomspace(atomspace,base_datasets_dir)
    total_patients = preprocess(atomspace)
    apply_subset_rule1(atomspace)
    apply_subset_rule2(atomspace)
    calculate_truth_values(atomspace, total_patients)
    remove_processed_subsets(atomspace)
    generate_attraction_links(atomspace)
    export_all_atoms(atomspace, base_results_dir)
    return atomspace

def parse_args():
    parser = argparse.ArgumentParser(description='Generate embedding vector for patients')
    parser.add_argument('--datapath', type=str, default='',
                        help='path to the source data')
    parser.add_argument('--outputpath', type=str, default='',
                        help='path to store the output files')
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

  if not os.path.exists(base_results_dir):
    os.makedirs(base_results_dir) 
  kb_as = generate_atoms(base_results_dir, base_datasets_dir)
  generate_embeddings("FMBPV",base_results_dir, kb_atomspace=kb_as)