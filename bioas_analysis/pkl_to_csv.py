import pandas as pd
import sys
from ast import literal_eval

def pkl_to_csv(pklfile):
    data = pd.read_pickle(pklfile, compression=None)

    df = pd.DataFrame([], columns=["patient_ID", "vector"])

    val = [list(i) for i in data.values()]
    vector_length = set([len(i) for i in val])
    df["patient_ID"] = data.keys()
    df["vector"] = val

    csv_file = "{}_expanded.csv".format(pklfile.replace(".pkl", ""))
    df = splitDataFrameListtocol(df, "vector")
    df.to_csv(csv_file, sep="\t", index=False)
    print("Vector length: {}, number of patients: {}".format(vector_length, df.shape[0]))

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

if __name__ == "__main__":
    pklfile = sys.argv[1]
    pkl_to_csv(pklfile)