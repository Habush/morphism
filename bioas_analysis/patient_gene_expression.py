__author__ = 'Abdulrahman Semrie<hsamireh@gmail.com>'

import pandas as pd
import argparse
import os
from atomwrappers import *
from zipfile import ZipFile
from io import BytesIO

def create_gene_expr_ln(patient, gene, value, fp, overexpr):
    if not pd.isna(value):
        if overexpr:
            gene_exec_ln = CLazyExecutionOutputLink(CSchemaNode("make-overexpression-schema-for-gene"), CGeneNode(gene))
            exec_ln = CExecutionLink(gene_exec_ln, patient,
                                     CNumberNode(str(value)))
        else:
            gene_exec_ln = CLazyExecutionOutputLink(CSchemaNode("make-underexpression-schema-for-gene"), CGeneNode(gene))
            exec_ln = CExecutionLink(gene_exec_ln, patient,
                                     CNumberNode(str(value)))

        fp.write(exec_ln.recursive_print() + "\n")


def expr_to_int(gene, val, expr_dict, overexpr):
    if overexpr:
        if val <= expr_dict[gene]: #For the overexpression,  we consider underexpression to the value 0
            return 0
        else:
            return val - expr_dict[gene]

    else:
        if val >= expr_dict[gene]: #For underexpression, we consider overexpression to have value 0
            return 0
        else:
            return expr_dict[gene] - val

def create_quantitative_predicate_ln(gene, overexpr, fp):
    if overexpr:
        quant_ln = CQuantitativePredicateLink(
                CLazyExecutionOutputLink(CSchemaNode("make-overexpression-schema-for-gene"), CGeneNode(gene)),
                CLazyExecutionOutputLink(CSchemaNode("make-overexpression-predicate-for-gene"), CGeneNode(gene)))
    else:
        quant_ln = CQuantitativePredicateLink(
            CLazyExecutionOutputLink(CSchemaNode("make-underexpression-schema-for-gene"), CGeneNode(gene)),
            CLazyExecutionOutputLink(CSchemaNode("make-underexpression-predicate-for-gene"), CGeneNode(gene)))

    fp.write(quant_ln.recursive_print() + "\n")

def create_schemanode(gene, patientnode, overexpr, fp, stv):
    if overexpr:
        lazy_ln = CLazyExecutionOutputLink(CSchemaNode("make-overexpression-predicate-for-gene"), CGeneNode(gene))
    else:
        lazy_ln = CLazyExecutionOutputLink(CSchemaNode("make-underexpression-predicate-for-gene"), CGeneNode(gene))
    eval_ln = CEvaluationLink(lazy_ln,patientnode, stv="(stv {} 1)".format(stv))
    fp.write(eval_ln.recursive_print() + "\n")

def create_member_ln(gene, fp):
    member_ln = CMemberLink(CGeneNode(gene), CConceptNode("profiled-genes"))
    fp.write(member_ln.recursive_print() + "\n")

def import_gene_expr(patient_df, overexpr_fp, underexpr_fp, genes_list):
    df = patient_df.dropna(axis=1, how='all')
    mean_dict = {}
    cols = df.columns
    # genes_list = [get_current_symbol(g) for g in genes_list]
    genes_list_ = [g for g in genes_list if g in cols]
    if len(genes_list_) < len(genes_list):
        dif = [i for i in genes_list if not i in genes_list_]
        print("No gene expression for the following genes: {}".format(dif))
    for col in df.loc[:, genes_list_]:
        if col == "patient_ID": continue
        mean_dict[str(col)] = df[col].median()

    for i in range(df.shape[0]):
        patient = CPatientNode(str(int(df.iloc[i]["patient_ID"])))
        for k in mean_dict:
            create_member_ln(k, overexpr_fp)
            create_member_ln(k, underexpr_fp)
            create_quantitative_predicate_ln(k, True, overexpr_fp)
            create_quantitative_predicate_ln(k, False, underexpr_fp)
            create_gene_expr_ln(patient, k, expr_to_int(k, df.iloc[i][k], mean_dict, True), overexpr_fp, True)
            create_gene_expr_ln(patient, k, expr_to_int(k, df.iloc[i][k], mean_dict, False), underexpr_fp, False)

def absgene_expr_atomese(patient_df, output_file, genes_list):
    cols = patient_df.columns
    genes_list_ = [g for g in genes_list if g in cols]
    if len(genes_list_) < len(genes_list):
        dif = [i for i in genes_list if not i in genes_list_]
        print("No gene expression for the following genes: {}".format(dif))
    df = patient_df[["patient_ID"] + genes_list_]
    df = df.dropna(axis=1, how='all')
    for col in df.columns:
        if col == "patient_ID": continue
        max_val = df[col].max()
        df[col] = df[col] / max_val
        create_member_ln(col, output_file)

    for i in range(df.shape[0]):
        patient = CPatientNode(str(int(df.iloc[i]["patient_ID"])))
        for gene in genes_list_:
            stv = "(stv {} 1.0)".format(df.iloc[i][gene])
            output_file.write(CEvaluationLink(CPredicateNode(gene),patient, stv=stv).recursive_print()+ "\n")

def qnormalized_atomse(df, overexpr_fp, underexpr_fp, genes):
    genes = set([g.replace("_underexpr","").replace("_overexpr","") for g in genes])
    for g in genes:
        create_member_ln(g, overexpr_fp)
        create_member_ln(g, underexpr_fp)
    for i in range(df.shape[0]):
        patient = CPatientNode(str(int(df.iloc[i]["patient_ID"])))
        for g in genes:
            try:
                underexpr_stv = df.iloc[i]["{}_underexpr".format(g)]
                create_schemanode(g, patient, False, underexpr_fp, underexpr_stv)
            except:
                pass
            try:    
                overexpr_stv = df.iloc[i]["{}_overexpr".format(g)]
                create_schemanode(g, patient, True,overexpr_fp, overexpr_stv)
            except:
                pass

def parse_args():
    parser = argparse.ArgumentParser(description="convert clinical trial patient data to atomese")
    parser.add_argument("--table", type=str, default='',
                        help="Path to clinical trial info table in csv format")
    parser.add_argument("--path", type=str, default='',
                        help="Path to save the output atomese file")
    parser.add_argument("--genes", type=str, default='',
                        help="Text file with list of selected genes")
    parser.add_argument("--pid", type=str, default='',
                        help="Text file with list of patient ID's")
    parser.add_argument("--relative", type=bool, default=False,
                        help="The gene expressions need to be normalized as under and over expression")
    parser.add_argument("--qnorm", type=bool, default=False,
                        help="The input dataset is quantile normalized")
    parser.add_argument("--prefix", type=str, default="",
                        help="The tag for the output file")
    parser.add_argument("--fqnorm", type=bool, default=False,
                        help="If this is true, all Gene expression values assumed to be over expressed")
    return parser.parse_args()

def download_data():
    merged_file_ln = "https://snet-bio-data.s3-us-west-2.amazonaws.com/example15bmc/merged-combat15.csv.xz"
    return build_request(merged_file_ln).read()

def main():
    print("Importing data")
    args = parse_args()
    if args.table:
        patient_df = pd.read_csv(args.table)
        if "True" in patient_df.columns:
            patient_df.rename({'True': 'patient_ID'}, axis=1, inplace=True)
    else:
        merged_file = download_data()
        patient_df = pd.read_csv(merged_file, compression="xz")

    if args.genes:
        genes_list = open(args.genes,"r").read().splitlines()
        patient_df = patient_df[["patient_ID"]+genes_list]
    else:
        genes_list = ["ERBB2", "GATA3", "TP53", "MKI67", "PIK3C3", "FOXA1", "BRCA2", "BRCA1"]

    if args.pid:
        pid = open(args.pid,"r").read().splitlines()
        patient_df = patient_df[patient_df["patient_ID"].isin(pid)]
    else:
        pid = False

    if args.path:
        save_path = os.path.abspath(args.path)
    else:
        save_path = os.getcwd()
    if args.prefix:
        prefix = args.prefix
    else:
        prefix = "all"
       
    if args.relative or args.qnorm or args.fqnorm:
        overexpr_file = os.path.join(save_path, "patient_gene_over_expr.scm")
        underexpr_file = os.path.join(save_path, "patient_gene_under_expr.scm")

        with open(overexpr_file, "w") as f1:
            with open(underexpr_file, "w") as f2:
                if args.qnorm:
                    print("relative median norm data")
                    qnormalized_atomse(patient_df, f1, f2, genes_list)
                elif args.fqnorm:
                    patient_df.columns = ["{}_overexpr".format(i) for i in patient_df.columns]
                    patient_df.rename({'patient_ID_overexpr': 'patient_ID'}, axis=1, inplace=True)
                    qnormalized_atomse(patient_df, f1, f2, genes_list)
                else:
                    print("Atomspace based qnorm will be applied later")
                    import_gene_expr(patient_df, f1, f2, genes_list)
    else:
        expr_file = os.path.join(save_path, "patient_gene_expr_absolute_{}.scm".format(prefix))
        with open(expr_file, "w") as f1:
            absgene_expr_atomese(patient_df, f1, genes_list)


if __name__ == "__main__":
    main()
    print("Done")