import os
import opencog.bioscience
from opencog.atomspace import AtomSpace, types
from opencog.execute import execute_atom
from opencog.logger import log
from opencog.scheme import scheme_eval
from opencog.type_constructors import *
from opencog.utilities import initialize_opencog
import pandas as pd 

df = pd.read_csv("cancer_data/bcTabs.csv", sep=",")

state_var = df[(df["varType"] == "state") & (df["subType"] == "preTx")]["variableName"].values

atomspace = AtomSpace()
initialize_opencog(atomspace)

scheme_eval(atomspace,""" 
(use-modules (opencog) (opencog bioscience))
(primitive-load "cancer_data/patient_data_all.scm")""")

evals = atomspace.get_atoms_by_type(types.EvaluationLink)

non_state_pred = []
for e in evals:
    pred = e.out[0].name
    if not pred in state_var:
        if not pred in non_state_pred:
            non_state_pred.append(pred)
        scheme_eval(atomspace, """(cog-delete {})""".format(e))

scheme_eval(atomspace, " ".join([
"(define (write-atoms-to-file file atoms)",
    "(define fp (open-output-file file))",
    "(for-each",
    "(lambda (x) (display x fp))",
    "atoms)",
    "(close-port fp))"]))

def write_atoms_to_file(filename, atom_list_str, atomspace):
    if not os.path.exists(filename):
        open(filename, "w").close()
        scheme_eval(atomspace, " ".join([
            "(write-atoms-to-file",
            "\"" + filename + "\"",
            atom_list_str + ")"]))

write_atoms_to_file("patient_data_state&preTX.scm", "(cog-get-atoms 'EvaluationLink)", atomspace)

print(non_state_pred)
