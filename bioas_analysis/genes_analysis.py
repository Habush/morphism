import os
from utils import *

log.set_level("ERROR")

def load_atomspace(datapath, atomspace):
  print("--- Load preprocessed data to the AtomSpace")
  scheme_eval(atomspace, "(load-file \"{}/{}\")".format(datapath, "member-links.scm"))  
  scheme_eval(atomspace, "(load-file \"{}/{}\")".format(datapath, "inheritance-links.scm"))
  # subset links are between GO catagories and we dont need them here 
  # scheme_eval(atomspace, "(load-file \"{}/{}\")".format(datapath, "subset-links.scm")) 

def generate_subsets(atomspace):
  print("--- Infer subsets between gene singletons and concepts")
  scheme_eval(atomspace, "(pln-load 'empty)")
  scheme_eval(atomspace, "(pln-load-from-path \"rules/subset-genes-rule.scm\")")
  scheme_eval(atomspace, "(pln-add-rule \"subset-genes-rule\")")
  scheme_eval(atomspace, " ".join([
    "(pln-bc (SubsetLink (Set (Variable \"$X\")) (Variable \"$Y\"))",
    "#:vardecl",
        "(VariableSet",
        "(TypedVariable (Variable \"$X\") (Type \"GeneNode\"))",
        "(TypedVariable (Variable \"$Y\") (TypeInh \"ConceptNode\")))",
    "#:maximum-iterations 12",
    "#:complexity-penalty 10)"]))

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
  scheme_eval(atomspace, "(pln-load 'empty)")
  scheme_eval(atomspace, "(pln-load-from-path \"rules/subset_negation_rule.scm\")")
  # (Subset (Set A) B) |- (Subset (Not (Set A) B))
  scheme_eval(atomspace, "(pln-add-rule \"subset-condition-negation-genes-rule\")")
  scheme_eval(atomspace, " ".join(["(pln-bc",
                  "(Subset (Not (Set (Variable \"$X\"))) (Variable \"$Y\"))",
                  "#:vardecl",
                    "(VariableSet",
                      "(TypedVariable (Variable \"$X\") (Type \"GeneNode\"))",
                      "(TypedVariable (Variable \"$Y\") (TypeInh \"ConceptNode\")))",
                  "#:maximum-iterations 12",
                  "#:complexity-penalty 10)"]))

def infer_attractions(atomspace):
  print("--- Inferring AttractionLinks")
  scheme_eval(atomspace, "(pln-load 'empty)")
  scheme_eval(atomspace, "(pln-load-from-path \"rules/subset_attraction_rule.scm\")")
  # (Subset (Set A) B)) (Subset (Not (Set A) B)) |- (Attraction (Set A) B)
  scheme_eval(atomspace, "(pln-add-rule \"subset-attraction-genes-rule\")")
  scheme_eval(atomspace, " ".join(["(pln-bc",
                  "(Attraction (Set (Variable \"$X\")) (Variable \"$Y\"))",
                  "#:vardecl",
                    "(VariableSet",
                      "(TypedVariable (Variable \"$X\") (Type \"GeneNode\"))",
                      "(TypedVariable (Variable \"$Y\") (TypeInh \"ConceptNode\")))",
                  "#:maximum-iterations 12",
                  "#:complexity-penalty 10)"]))

def export_result(datapath, atomspace):
  attraction_links_scm = os.path.join(datapath, "attraction-links-genes.scm")
  print("--- Exporting attraction links to files")
  scheme_eval(atomspace, " ".join([
  "(define (write-atoms-to-file file atoms)",
    "(define fp (open-output-file file))",
    "(for-each",
      "(lambda (x) (display x fp))",
      "atoms)",
    "(close-port fp))"]))
  write_atoms_to_file(attraction_links_scm, "(cog-get-atoms 'AttractionLink)", atomspace)

def generate_attractionlinks_genes(datapath):
  ### Initialize the AtomSpace ###
  atomspace = AtomSpace()
  initialize_opencog(atomspace)

  scheme_eval(atomspace, "(add-to-load-path \".\")")
  scheme_eval(atomspace, "(use-modules (opencog) (opencog bioscience) (opencog ure) (opencog pln) (opencog persist-file) (srfi srfi-1))")
  load_atomspace(datapath, atomspace)
  gene_singleton_tv(atomspace)
  generate_subsets(atomspace)
  infer_negation(atomspace)
  infer_attractions(atomspace)
  export_result(datapath, atomspace)
