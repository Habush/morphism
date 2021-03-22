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

zerovector = []

def generate_embeddings(data_dir, node_type, kb_atomspace=False, p=1, T=2,normalization=False,
                        test_data=False, vector_space=False):
  if kb_atomspace:
    atomspace = kb_atomspace
  else:
    print("--- Loading the AtomSpace")
    atomspace = AtomSpace()
    initialize_opencog(atomspace)
    scheme_eval(atomspace, "(use-modules (opencog persist-file))")
    for i in os.listdir(data_dir):
        if i.endswith(".scm"):
            scheme_eval(atomspace,"(load-file \"{}/{}\")".format(data_dir, i))

  property_vectors = build_property_vectors(atomspace, data_dir, node_type, p, T, normalization=normalization, 
                                            test_data=test_data, vector_space=vector_space)
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
    embedding_vector = do_kpca(property_vectors, test_data=test_data, vector_space=vector_space)
  else:    
    embedding_vector = do_kpca(property_vectors)

  result = export_property_vectors(data_dir, embedding_vector,p,T, normalization)
  return result

def build_property_vectors(atomspace, data_dir, node_type, p, T, normalization=False,test_data=False, vector_space=False):
  print("--- Building property vectors")
  ppty = set([i.out[1] for i in atomspace.get_atoms_by_type(types.AttractionLink) if i.out[1].type_name != "VariableNode"])
  nodes = set([i.out[0] for i in atomspace.get_atoms_by_type(types.AttractionLink)])
  ppty_mapping = dict(zip(range(len(ppty)), [str(i) for i in ppty]))

  property_df = pd.DataFrame([], columns=["patient_ID"] + list(ppty_mapping.keys()))
  print("Number of {}: {}".format(node_type, len(atomspace.get_atoms_by_type(getattr(types, node_type)))))
  print("Number of properties: {}".format(len(ppty)))
  for node in nodes:
    if node.type_name == node_type or node.out[0].type_name == node_type:
      p_vec = []
      for pt in ppty:
          if pt.atomspace.is_link_in_atomspace(types.AttractionLink, [node, pt]):
              attraction = AttractionLink(node, pt)
              attraction_tv = attraction.tv
              # weighted product, to give more weight to the strength s^p * c^(T-p)
              # with default p = 1 and T = 2
              wp = attraction.tv.mean**p * attraction.tv.confidence**(T-p)
              p_vec.append(wp)
          else:
              p_vec.append(0.0)

      if sum(p_vec) != 0:
        node_name = node.out[0].name if node.type_name == "SetLink" else node.name 
        property_df.loc[len(property_df),:] = [node_name] + p_vec
      else:
        zerovector.append(str(node))
        continue

  # normalize the vector
  if normalization:
    if test_data:
        property_df = normalize_df(property_df, normalization,test_data=True, train_data=vector_space)
    else:
        property_df = normalize_df(property_df, normalization)

  # Dump the property vector to CSV (before applying kpca)
  property_vector_pickle_beforekpca = os.path.join(data_dir, 
            "property_vector_beforekpca_p={},T={}_{}_{}.csv".format(p, T, normalization if normalization else "notnormalized",str(date.today())))
  property_df.to_csv(property_vector_pickle_beforekpca, sep="\t", index=False)
  filehandler = open(os.path.join(data_dir, "ppty_mapping_{}".format(normalization if normalization else "notnormalized")), 'wb')
  pickle.dump(ppty_mapping, filehandler)
  return property_df

def normalize_df(property_df, normalization, test_data=False, train_data=False):
  cols = property_df.columns.drop("patient_ID")
  if test_data:
    train_data = train_data.drop("patient_ID")
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
    qt = QuantileTransformer(n_quantiles=5, random_state=0, output_distribution="normal")
    if test_data:
      qt.fit(train_data)
      property_df[cols] = qt.transform(property_df[cols])
    else:
      property_df[cols] = qt.fit_transform(property_df[cols])
  return property_df

def do_kpca(property_vectors_df, test_data=False, vector_space=False, n_components=150, tan_dist=False):
  def kernel_func(X, Y):
    X = sparse.vstack(X.values())
    Y = sparse.vstack(Y.values())
    len_x = X.shape[0]
    len_y = Y.shape[0]
    dist_array = numpy.random.random((len_x, len_y)).astype(numpy.float16) * 0 - 1.0
    i = 0
    for a in X:
      a = a.toarray().flatten()
      j = 0
      for b in Y:
        b = b.toarray().flatten()
        if tan_dist:
          dist = tanimoto_kernel(a, b)
        else:
          dist = tanimoto(a, b)
        dist_array[i, j] = dist
        j += 1
      i += 1
    return dist_array

  print("--- Doing KPCA")
  kpca = KernelPCA(kernel = "precomputed", n_components=n_components)
  property_vectors_dict = compress_df(property_vectors_df)
  if test_data:
    print("Train data:{}, Test data:{}".format(vector_space.shape, property_vectors_df.shape))
    train_data = compress_df(vector_space)
    train_space = kernel_func(train_data, train_data)
    # Fit with train data and do just transform on test
    kpca.fit(train_space)
    test_space = kernel_func(property_vectors_dict, train_data)
    print("--- Test data : {}".format(test_space.shape))
    result_kpca = kpca.transform(test_space)

  else:
    result_kpca = kpca.fit_transform(kernel_func(property_vectors_dict, property_vectors_dict))
 
  for k, kpca_v in zip(property_vectors_dict.keys(), result_kpca):
    property_vectors_dict[k] = kpca_v
  return dict_to_df(property_vectors_dict)

def dict_to_df(res_dict):
  df = pd.DataFrame(res_dict.values())
  df["patient_ID"] = res_dict.keys()
  return df

def compress_df(df):
  # Compress the dataframe (sparse matrix) to dict
  vectors_dict = {}
  for i,j in enumerate(df["patient_ID"]):
    vectors_dict[j] = sparse.csr_matrix(df.loc[i,df.columns.drop("patient_ID")].values.tolist())
  return vectors_dict

def export_property_vectors(data_dir, property_vectors,p,T, normalization):
  if not normalization:
    normalization = "notnormalized"
  output_file = os.path.join(data_dir, 
          "property_vector_afterkpca_p={},T={}_{}_{}.csv".format(p,T,normalization,str(date.today())))

  print("--- Exporting property vectors to \"{}\"".format(output_file))

  property_vectors.to_csv(output_file, sep="\t", index=False)

  if len(zerovector) > 0:
    print(len(zerovector))
    with open(os.path.join(data_dir , "zerovector.txt"), "w") as z:
      z.write("\n".join(zerovector))
  return property_vectors

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
