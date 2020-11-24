import os
# from preprocess import *
from utils import *

def load_atomspace(datapath, atomspace):
  print("--- Load preprocessed data to the AtomSpace")
  scheme_eval(atomspace, "(load-file \"{}/{}\")".format(datapath, "member-links.scm"))  
  scheme_eval(atomspace, "(load-file \"{}/{}\")".format(datapath, "inheritance-links.scm"))  
  scheme_eval(atomspace, "(load-file \"{}/{}\")".format(datapath, "subset-links.scm")) 

def generate_direct_subsets(atomspace):
  print("--- Infer direct subsets between concepts")
  scheme_eval(atomspace, "(pln-load 'empty)")
  scheme_eval(atomspace, "(pln-load-from-path \"/home/hedra/OpenCOg/pln/opencog/pln/rules/extensional/subset-direct-introduction.scm\")")
  scheme_eval(atomspace, "(pln-add-rule \"subset-direct-introduction-rule\")")
  scheme_eval(atomspace, " ".join([
    "(pln-bc (SubsetLink (Variable \"$X\") (Variable \"$Y\"))",
    "#:vardecl",
        "(VariableSet",
        "(TypedVariable (Variable \"$X\") (TypeInh \"ConceptNode\"))",
        "(TypedVariable (Variable \"$Y\") (TypeInh \"ConceptNode\")))",
    "#:maximum-iterations 12",
    "#:complexity-penalty 10)"]))

def infer_attractions(atomspace):
  print("--- Inferring AttractionLinks")
  scheme_eval(atomspace, "(pln-load 'empty)")
  scheme_eval(atomspace, "(pln-load-from-path \"/home/hedra/OpenCOg/pln/opencog/pln/rules/intensional/attraction-introduction.scm\")")
  scheme_eval(atomspace, "(pln-load-from-path \"/home/hedra/OpenCOg/pln/opencog/pln/rules/term/condition-negation.scm\")")
  # (Subset A B) |- (Subset (Not A) B)
  scheme_eval(atomspace, "(pln-add-rule \"subset-condition-negation-rule\")")
  # (Subset A B) (Subset (Not A) B) |- (Attraction A B)
  scheme_eval(atomspace, "(pln-add-rule \"subset-attraction-introduction-rule\")")
  scheme_eval(atomspace, " ".join(["(pln-bc",
                  "(Attraction (Variable \"$X\") (Variable \"$Y\"))",
                  "#:vardecl",
                    "(VariableSet",
                      "(TypedVariable (Variable \"$X\") (TypeInh \"ConceptNode\"))",
                      "(TypedVariable (Variable \"$Y\") (TypeInh \"ConceptNode\")))",
                  "#:maximum-iterations 12",
                  "#:complexity-penalty 10)"]))

def export_result(datapath, atomspace):
  attraction_links_scm = os.path.join(datapath, "attraction-links.scm")
  print("--- Exporting attraction links to files")
  write_atoms_to_file(attraction_links_scm, "(cog-get-atoms 'AttractionLink)", atomspace)

def generate_attractionlinks_concepts(datapath):
  ### Initialize the AtomSpace ###
  atomspace = AtomSpace()
  initialize_opencog(atomspace)

  scheme_eval(atomspace, "(add-to-load-path \".\")")
  scheme_eval(atomspace, "(use-modules (opencog) (opencog bioscience) (opencog ure) (opencog pln) (opencog persist-file) (srfi srfi-1))")
  scheme_eval(atomspace, " ".join([
    "(define (write-atoms-to-file file atoms)",
      "(define fp (open-output-file file))",
      "(for-each",
        "(lambda (x) (display x fp))",
        "atoms)",
      "(close-port fp))"]))
  load_atomspace(datapath, atomspace)
  generate_direct_subsets(atomspace)
  infer_attractions(atomspace)
  export_result(datapath, atomspace)
