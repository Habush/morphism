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

def dict_to_csv(pklfile, file_name=False):
    if file_name:
        data = pd.read_pickle(pklfile, compression=None)
    else:
        data = pklfile
    df = pd.DataFrame([], columns=["patient_ID", "vector"])
    val = [list(i) for i in data.values()]
    df["patient_ID"] = data.keys()
    df["vector"] = val
    df = splitDataFrameListtocol(df, "vector")
    if file_name:
        csv_file = "{}_expanded.csv".format(pklfile.replace(".pkl", ""))
        df.to_csv(csv_file, sep="\t", index=False)
    print("Output vector shape: {}".format(df.shape))
    return df

def splitDataFrameListtocol(df,target_column):
    ''' df = dataframe to split,
    target_column = the column containing the values to split
    separator = the symbol used to perform the split
    returns: a dataframe with each entry for the target column separated and moved into a new column. 
    '''
    # all columns except `target_column`
    others = df.columns.difference([target_column])
    new_df = df[others]

    data_expanded = [str(i).replace("]","").replace("[","").split(",") for i in df[target_column].tolist()]
    df2 = pd.DataFrame(data_expanded, columns=range(len(data_expanded[0])))

    return pd.concat([new_df, df2], axis=1)

def quantile_normalize(df):
    patient_id = df["patient_ID"]
    df = df.drop("patient_ID", axis=1)
    """
    input: dataframe with numerical columns
    output: dataframe with quantile normalized values
    """
    df_sorted = pd.DataFrame(np.sort(df.values,
                                     axis=0), 
                             index=df.index, 
                             columns=df.columns)
    df_mean = df_sorted.mean(axis=1)
    df_mean.index = np.arange(1, len(df_mean) + 1)
    df_qn =df.rank(method="min").stack().astype(int).map(df_mean).unstack()
    df_qn["patient_ID"] = patient_id
    return(df_qn)