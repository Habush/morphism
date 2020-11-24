import os
import opencog.bioscience
from opencog.atomspace import AtomSpace, types
from opencog.execute import execute_atom
from opencog.logger import log
from opencog.scheme import scheme_eval
from opencog.type_constructors import *
from opencog.utilities import initialize_opencog

log.set_level("ERROR")


## Helper Function and variables
def get_concepts(atomspace, lst):
    result = []
    for i in lst:
        result = result + atomspace.get_atoms_by_type(getattr(types, i))
    return result 

concepts_lst = ["CellularComponentNode","BiologicalProcessNode","MolecularFunctionNode","ReactomeNode","CellNode","UberonNode"]
genes_lst = ["GeneNode", "SetLink"]

def write_atoms_to_file(filename, atom_list_str, atomspace):
    if not os.path.exists(filename):
        open(filename, "w").close()
        scheme_eval(atomspace, " ".join([
            "(write-atoms-to-file",
            "\"" + filename + "\"",
            atom_list_str + ")"]))