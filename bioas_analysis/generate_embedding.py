import numpy
from scipy import sparse
import os
import json
import pickle
import random
from gensim.models import Word2Vec
from matplotlib import pyplot
from scipy.spatial import distance
from scipy.stats import kendalltau, pearsonr, spearmanr
from sklearn.decomposition import PCA, KernelPCA
from sklearn.preprocessing import StandardScaler, MinMaxScaler, QuantileTransformer
from datetime import date
from utils import *
import pandas as pd
import re

zerovector = []
nodetypes = ["PredicateNode","BiologicalProcessNode", "ReactomeNode", "GeneNode","MolecularFunctionNode","CellularComponentNode","PharmGkbNode"]

def generate_embeddings(data_dir, node_type, kb_atomspace=False, ppty_vector=False, p=1, T=2,normalization=False,
                        test_data=False, vector_space=False, abs_ge=False):
  if kb_atomspace:
    atomspace = kb_atomspace
    property_vectors = build_property_vectors(atomspace, node_type, p, T, abs_ge)
  else:
    if ppty_vector:
      property_vectors = pd.read_csv(ppty_vector, sep="\t")
      print("Property vector shape {}".format(property_vectors.shape))
    else:
      print("--- Loading the AtomSpace")
      atomspace = AtomSpace()
      initialize_opencog(atomspace)
      scheme_eval(atomspace, "(use-modules (opencog persist-file))")
      scheme_eval(atomspace,"(load-file \"{}/{}\")".format(data_dir, "AttractionLinks-results.scm.bc"))
      property_vectors = build_property_vectors(atomspace, node_type, p, T, abs_ge)

  # normalize the vector
  if normalization:
    if test_data:
        property_vectors = normalize_df(property_vectors, normalization,test_data=True, train_data=vector_space)
    else:
        property_vectors = normalize_df(property_vectors, normalization)

  if not ppty_vector:
  # Dump the property vector to CSV (before applying kpca)
    property_vector_beforekpca = os.path.join(data_dir, 
              "property_vector_beforekpca_p={},T={}_{}_{}.csv".format(p, T, normalization if normalization else "notnormalized",str(date.today())))
    property_vectors.to_csv(property_vector_beforekpca, sep="\t", index=False)
  
  node_ID = "patient_ID" if node_type == "PatientNode" else "node_ID"
  if test_data:
    # Patients in the test or training might not have all attributes, fill missing ones with 0
    if len(property_vectors.columns) < len(vector_space.columns):
      rows = property_vectors.shape[0]
      col_diff = len(vector_space.columns) - len(property_vectors.columns)
      for i in range(col_diff):
        property_vectors["missign{}".format(i)] = [0]*rows
    elif len(property_vectors.columns) > len(vector_space.columns):
      rows = vector_space.shape[0]
      col_diff = len(property_vectors.columns) - len(vector_space.columns)
      for i in range(col_diff):
        vector_space["missign{}".format(i)] = [0]*rows
    embedding_vector = do_kpca(property_vectors.drop(node_ID, axis=1), test_data=test_data, vector_space=vector_space.drop(node_ID, axis=1))
  else:    
    embedding_vector = do_kpca(property_vectors.drop(node_ID, axis=1))
  embedding_vector[node_ID] = property_vectors[node_ID]
  result = export_property_vectors(data_dir, embedding_vector,p,T, normalization)
  return result

def get_node(satLink, abs_ge=False):
    result = satLink
    nodes = re.findall('\(.*?\)',satLink)
    if nodes:
      for i in nodes:
        if any(a in i for a in nodetypes):
            result = i
            break
    postfix = ""
    if "overexpression" in satLink:
        postfix = "_overexp"
    elif "underexpression" in satLink:
        postfix = "_underexp"
    if abs_ge:
        postfix = ""
    result = "{}{}".format(re.findall('"([^"]*)"', result)[0], postfix)
    return result

def build_property_vectors(atomspace, node_type, p, T, abs_ge):
  print("--- Building property vectors")
  ppty = set([i.out[1] for i in atomspace.get_atoms_by_type(types.AttractionLink) if i.out[1].type_name != "VariableNode"])
  nodes = set([i.out[0] for i in atomspace.get_atoms_by_type(types.AttractionLink)])
  print("Number of {}: {}".format(node_type, len(nodes)))
  print("Number of properties: {}".format(len(ppty)))
  property_df = pd.DataFrame(numpy.zeros((len(nodes), len(ppty))), columns=ppty)
  node_ID = "patient_ID" if node_type == "PatientNode" else "node_ID"
  property_df[node_ID] = nodes
  property_df = property_df.set_index(node_ID)
  for att_link in atomspace.get_atoms_by_type(types.AttractionLink):
    node = att_link.out[0]
    pt = att_link.out[1]
    attraction_tv = att_link.tv
    # weighted product, to give more weight to the strength s^p * c^(T-p)
    # with default p = 1 and T = 2
    wp = attraction_tv.mean**p * attraction_tv.confidence**(T-p)
    property_df.loc[node][pt] = wp

  property_df.columns = [get_node(str(i), abs_ge=abs_ge) for i in property_df.columns]
  property_df = property_df.loc[~(property_df==0).all(axis=1)] # remove zero vector
  property_df.reset_index(inplace = True)
  property_df[node_ID] = [get_node(str(i), abs_ge=abs_ge) for i in property_df[node_ID]]
  return property_df

def normalize_df(property_df, normalization, test_data=False, train_data=False):
  cols = property_df.columns.drop("patient_ID")
  if test_data:
    train_data = train_data.drop("patient_ID", axis=1)
  if normalization == "scaling":
    scaler = MinMaxScaler()
    if test_data:
      scaler.fit(train_data)
      property_df[cols] = scaler.transform(property_df[cols])
    else:
      property_df[cols] = scaler.fit_transform(property_df[cols])
  elif normalization == "standard":
    scaler = StandardScaler()
    if test_data:
      scaler.fit(train_data)
      property_df[cols] = scaler.transform(property_df[cols])
    else:
      property_df[cols] = scaler.fit_transform(property_df[cols])
  else:
    qt = QuantileTransformer(n_quantiles=5, random_state=0)
    if test_data:
      qt.fit(train_data)
      property_df[cols] = qt.transform(property_df[cols])
    else:
      property_df[cols] = qt.fit_transform(property_df[cols])
  return property_df

def tanimoto(v1, v2):
  v1_v2 = numpy.dot(v1, v2)
  v1_sq = numpy.sum(numpy.square(v1))
  v2_sq = numpy.sum(numpy.square(v2))
  return v1_v2 / (v1_sq + v2_sq - v1_v2)

def tanimoto_kernel(v1, v2):
    v1_norm = numpy.linalg.norm(v1, ord=1)
    v2_norm = numpy.linalg.norm(v2, ord=1)
    v1_v2_norm = numpy.linalg.norm((v1 - v2), ord=1)
    f = (v1_norm + v2_norm - v1_v2_norm)
    g = (v1_norm + v2_norm + v1_v2_norm)
    result = f / g
    return result

def intensional_similarity(v1, v2):
  result = sum(min(v1, v2))/sum(max(v1, v2))
  return result

def do_kpca(property_vectors_df, test_data=False, vector_space=False, n_components=150, kernel=tanimoto):
  print("--- Doing KPCA")
  kpca = KernelPCA(kernel=kernel, n_components=n_components, n_jobs=-1)
  if test_data:
    print("Train data:{}, Test data:{}".format(vector_space.shape, property_vectors_df.shape))
    # Fit with train data and do just transform on test
    kpca.fit(vector_space)
    result_kpca = kpca.transform(property_vectors_df)
  else:
    result_kpca = kpca.fit_transform(property_vectors_df)
 
  result_df = pd.DataFrame(result_kpca)
  return result_df

def export_property_vectors(data_dir, property_vectors,p,T, normalization):
  if not normalization:
    normalization = "notnormalized"
  output_file = os.path.join(data_dir, 
          "embedding_vector_p={},T={}_{}_{}.csv".format(p,T,normalization,str(date.today())))
  print("--- Exporting property vectors to \"{}\"".format(output_file))
  property_vectors.to_csv(output_file, sep="\t", index=False)
  if len(zerovector) > 0:
    print(len(zerovector))
    with open(os.path.join(data_dir , "zerovector.txt"), "w") as z:
      z.write("\n".join(zerovector))
  return property_vectors

def apply_kernel(X, Y, kernel_func=tanimoto):
    if X.shape[1] != Y.shape[1]: #we assume the two inputs to have the same dimensions
      raise ValueError("Different dimensions, dim(X): {0}, dim(Y): {1}".format(X.shape[1], Y.shape[1]))
    K = np.zeros((X.shape[0], Y.shape[0])) #kernel matrix
    i = 0
    for x in X:
        a = x.flatten()
        j = 0
        for y in Y:
            b = y.flatten()
            K[i, j] = kernel_func(a,b)
            j += 1        
        i += 1    
    return K