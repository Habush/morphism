import os
import opencog.bioscience
from opencog.atomspace import AtomSpace, types
from opencog.execute import execute_atom
from opencog.logger import log
from opencog.scheme import scheme_eval
from opencog.type_constructors import *
from opencog.utilities import initialize_opencog
import pandas as pd
import numpy as np
import sys

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

def filter_bp(atomspace):
    sub = atomspace.get_atoms_by_type(types.SubsetLink)

    # --- Example subset
    # (SubsetLink (stv 0.00619435 0.171843)
    # (AndLink
    #     (ConceptNode "profiled-genes")
    #     (BiologicalProcessNode "GO:0071695"))
    # (SatisfyingSetScopeLink
    #     (VariableNode "$G")
    #     (EvaluationLink
    #     (LazyExecutionOutputLink
    #         (SchemaNode "make-overexpression-predicate-for-gene")
    #         (VariableNode "$G"))
    #     (PatientNode "615289"))))
    # ----
    for s in sub:
        subelem = s.out[0]
        if subelem.type_name == "AndLink":
            elem1 = subelem.out[0]
            elem2 = subelem.out[1]
            if not str(elem1.type_name) != "BiologicalProcessNode" or str(elem2.type_name) != "BiologicalProcessNode":
                scheme_eval(atomspace, "(cog-delete! {})".format(s))    

    return atomspace
