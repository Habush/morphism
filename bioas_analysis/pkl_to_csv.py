import pandas as pd
import sys

def pkl_to_csv(pklfile):
    data = pd.read_pickle(pklfile, compression=None)

    df = pd.DataFrame([], columns=["patient_id", "vector"])

    val = [list(i) for i in data.values()]
    vector_length = set([len(i) for i in val])
    df["patient_id"] = data.keys()
    df["vector"] = val

    csv_file = "{}.csv".format(pklfile.replace(".pkl", ""))
    df.to_csv(csv_file, sep="\t", index=False)
    print("Vector length: {}, number of patients: {}".format(vector_length, df.shape[0]))

if __name__ == "__main__":
    pklfile = sys.argv[1]
    pkl_to_csv(pklfile)