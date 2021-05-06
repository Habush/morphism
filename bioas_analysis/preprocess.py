import os
from datetime import date
from utils import *

base_datasets_dir = os.getcwd() + "/bioas/"
base_results_dir = os.getcwd() + "/results/bioas-{}/".format(str(date.today()))
if not os.path.exists(base_results_dir):
  os.makedirs(base_results_dir)
member_links_scm = base_results_dir + "member-links.scm"
inheritance_links_scm = base_results_dir + "inheritance-links.scm"
subset_links_scm = base_results_dir + "subset-links.scm"

def populate_atomspace(atomspace):
  print("--- Populating the AtomSpace")
#   for i in os.listdir(base_datasets_dir):
#     scheme_eval(atomspace, "(load-file \"{}/{}\")".format(base_datasets_dir, i))
  scheme_eval(atomspace, "(load-file \"sample.scm\")")   

def partof_to_inheritance(atomspace):
  print("--- Converting has_part to Inheritance")
  # (Evaluation Predicate "has_part" List C1 C2) |- (Inheritance C2 C1)
  for e in [i for i in atomspace.get_atoms_by_type(types.EvaluationLink) if i.out[0].name == "has_part"]:
      lst = e.out[1].out
      scheme_eval(atomspace, "(InheritanceLink {} {})".format(lst[1], lst[0]))

def infer_subsets(atomspace):
  scheme_eval(atomspace, "(add-to-load-path \".\")")
  scheme_eval(atomspace, "(load-from-path \"rules/translation.scm\")")
  scheme_eval(atomspace, "(load-from-path \"rules/transitivity.scm\")")
  # (Inheritance C1 C2) |- (Subset C1 C2)
  print("--- Applying Translation Rule")
  scheme_eval(atomspace, "(inheritance->subset)")
  print("--- Applying Transitivity Rules")
  # (Subset C1 C2) (Subset C2 C3) |- (Subset C1 C3)
  scheme_eval(atomspace, "(gen-present-link-transitivity)")
  # (Member G C1) (Subset C1 C2) |- (Member G C2)
  scheme_eval(atomspace, "(gen-present-mixed-link-transitivity)")

def calculate_truth_values(atomspace):
  print("--- Calculating Truth Values")

  node_member_dict = {}
  def get_members(node):
    if node_member_dict.get(node):
      return node_member_dict[node]
    else:
      members = [x.out[0] for x in node.incoming if x.type == types.MemberLink and x.out[1] == node]
      node_member_dict[node] = members
      return members

  def get_confidence(count):
    return float(scheme_eval(atomspace, "(count->confidence " + str(count) + ")"))

  # MemberLinks (without pre-assigned TV) directly generated from the data, can be considered as true
  for m in atomspace.get_atoms_by_type(types.MemberLink):
    if m.tv.confidence == 0.0:
        m.tv = TruthValue(1, 1)

  # ConceptNode "A" (stv s c)
  # where:
  # s = |A| / |universe|
  # c = |universe|
  universe_size = len(atomspace.get_atoms_by_type(types.GeneNode))
  tv_confidence = get_confidence(universe_size)
  for c in get_concepts(atomspace, concepts_lst):
    member_size = len(get_members(c))
    tv_strength = member_size / universe_size
    c.tv = TruthValue(tv_strength, tv_confidence)

  # SubLinks (infered from inheritance link) are directly from the data, and is true by definition
  for s in atomspace.get_atoms_by_type(types.SubsetLink):
    s.tv = TruthValue(1, 1)

def export_all_atoms(atomspace):
  print("--- Exporting Atoms to files")
  write_atoms_to_file(member_links_scm, "(cog-get-atoms 'MemberLink)", atomspace)
  write_atoms_to_file(inheritance_links_scm, "(cog-get-atoms 'InheritanceLink)", atomspace)
  write_atoms_to_file(subset_links_scm, "(cog-get-atoms 'SubsetLink)", atomspace)

def generate_atoms():
    ### Initialize the AtomSpace ###
    atomspace = AtomSpace()
    initialize_opencog(atomspace)

    ### Guile setup ###
    scheme_eval(atomspace, "(add-to-load-path \".\")")
    scheme_eval(atomspace, "(use-modules (opencog) (opencog bioscience) (opencog ure) (opencog pln) (opencog logger) (opencog persist-file) (srfi srfi-1) (ice-9 threads))")
    scheme_eval(atomspace, " ".join([
    "(define (write-atoms-to-file file atoms)",
        "(define fp (open-output-file file))",
        "(for-each",
        "(lambda (x) (display x fp))",
        "atoms)",
        "(close-port fp))"]))

    populate_atomspace(atomspace)
    partof_to_inheritance(atomspace)
    infer_subsets(atomspace)
    calculate_truth_values(atomspace)
    export_all_atoms(atomspace)
    return base_results_dir

if __name__ == "__main__":
    output_path = generate_atoms()
    print("Output path {}".format(output_path))
