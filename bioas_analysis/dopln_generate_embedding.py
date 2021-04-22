from patient_gene_expression import *
from patients_analysis import *
from generate_embedding import build_property_vectors
from multiprocessing import Process
import shutil

def do_pln(dataset_path, overexpr, output_path):
    # install https://github.com/Habush/pln-brca-xp.git for this function to work
    atomspace = AtomSpace()
    initialize_opencog(atomspace)
    scheme_eval(atomspace, "(use-modules (pln-bio bio-utils) (pln-bio expr) (pln-bio concept-deduction))")
    if overexpr:
        print("--- Overexpressed genes PLN inference")
        scheme_eval(atomspace, """(run-concept-deduction-expr #t "kbs/plndata_dir" "{}")""".format(output_path))
        output_file = os.path.join(output_path, "../subset_over_expr.scm")
        concat_files(os.path.join(output_path, "batches_overexpr"), output_file)
    else:
        print("--- Underexpressed genes PLN inference")
        scheme_eval(atomspace, """(run-concept-deduction-expr #f "kbs/plndata_dir" "{}")""".format(output_path))
        output_file = os.path.join(output_path, "../subset_under_expr.scm")
        concat_files(os.path.join(output_path, "batches_underexpr"), output_file)

def concat_files(dataset_path, output_file):
    filenames = os.listdir(dataset_path)
    with open(output_file, 'w') as outfile:
        for fname in filenames:
            with open(os.path.join(dataset_path, fname)) as infile:
                for line in infile:
                    outfile.write(line)

def parse_arguments():
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
                        help="Tells if the gene expression values are normalized as under and over expression")
    parser.add_argument("--scaled", type=bool, default=False,
                        help="Tells if the gene expression values are scaled between 0 and 1 to be used as TV (for non-normalized data)")
    parser.add_argument("--infergo", type=bool, default=False,
                        help="If True, PLN inference to get GO/pathway for patients will be done")
    parser.add_argument("--plnonly", type=bool, default=False,
                        help="If True, only PLN inferred GO/pathway of patients will be used to generate the embedding")
    return parser.parse_args()

def create_dir(path):
    if not os.path.isdir(path):
        os.mkdir(path)

if __name__ == "__main__":
    args = parse_arguments()
    data = pd.read_csv(args.table)
    if args.pid:
        pid = open(args.pid,"r").read().splitlines()
        data = data[data["patient_ID"].isin(pid)]
    if args.genes:
        genes_list = open(args.genes,"r").read().splitlines()
    else:
        genes_list = data.columns.drop("patient_ID")
    if args.path:
        save_path = os.path.abspath(args.path)
    else:
        save_path = os.getcwd()
    if args.infergo:
        # pln_inference load-kbs function required the data to be listed under the same dir
        # replace "/home/hedra/pln-brca-xp/kbs" by where you installed pln-brca-xp 
        # copy the background knowlegde base from https://mozi.ai/datasets/kbs/
        plnoutput_path = os.path.join(save_path, "plnoutput_dir")
        plndata_path = os.path.join("/home/hedra/pln-brca-xp/kbs" ,"scm_data_dir")
        create_dir(plnoutput_path)
    else:
        plndata_path = os.path.join(save_path ,"scm_data_dir")
    create_dir(plndata_path)
    print("--- Convert gene expression values to Atomese format")
    # output files
    overexpr_file = os.path.join(plndata_path, "patient_gene_over_expr.scm")
    underexpr_file = os.path.join(plndata_path, "patient_gene_under_expr.scm")
    with open(overexpr_file, "w") as f1:
        with open(underexpr_file, "w") as f2:
            if args.relative:
                abs_ge = False
                qnormalized_atomse(data, f1, f2, genes_list)
            else:
                abs_ge = True
                if not args.scaled:
                    data = normalize_df(data, "scaling")
                data.columns = ["{}_overexpr".format(i) for i in data.columns]
                data.rename({'patient_ID_overexpr': 'patient_ID'}, axis=1, inplace=True)
                qnormalized_atomse(data, f1, f2, genes_list)
    if args.infergo:
        print("--- Doing PLN inference")
        if abs_ge:
            do_pln(plndata_path, True, plnoutput_path)
        else:   
            p1 = Process(target=do_pln, args=(plndata_path, True, plnoutput_path))
            p1.start()
            p2 = Process(target=do_pln, args=(plndata_path, False, plnoutput_path))
            p2.start()
            p1.join()
            p2.join()
        shutil.move(plndata_path, save_path)           
    print("--- Generate patient embeddings")
    if args.plnonly:
        kb_as = generate_atoms(save_path, save_path, False)
    else:
        data_src = os.path.join(save_path ,"scm_data_dir")
        for sub in os.listdir(save_path):
            if sub.endswith(".scm"):
                shutil.move(sub, data_src)
        kb_as = generate_atoms(save_path, data_src, False)
    export_all_atoms(kb_as, save_path)
    generate_embeddings(save_path,"PatientNode", kb_atomspace=kb_as, abs_ge=abs_ge)
    print("Done")