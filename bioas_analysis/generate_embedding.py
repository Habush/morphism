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
from datetime import date
from utils import *
import pandas as pd

zerovector = []

def generate_embeddings(embedding_method, data_dir, node_type, kb_atomspace=False, p=1.5, T=2,normalization=False):
    if embedding_method == "FMBPV":
      if kb_atomspace:
        atomspace = kb_atomspace
      else:
        print("--- Load the AtomSpace")
        atomspace = AtomSpace()
        initialize_opencog(atomspace)
        scheme_eval(atomspace, "(use-modules (opencog persist-file))")
        for i in os.listdir(data_dir):
            if i.endswith(".scm"):
                scheme_eval(atomspace,"(load-file \"{}/{}\")".format(data_dir, i))

      property_vectors = build_property_vectors(atomspace, data_dir, node_type, p, T, normalization=normalization)
      property_vectors = do_kpca(property_vectors)
      result = export_property_vectors(data_dir, property_vectors,p,T, normalization)
      return result

def build_property_vectors(atomspace, data_dir, node_type, p, T, normalization=False):
  print("--- Building property vectors")
  ppty = set([i.out[1] for i in atomspace.get_atoms_by_type(types.AttractionLink) if i.out[1].type_name != "VariableNode"])
  nodes = set([i.out[0] for i in atomspace.get_atoms_by_type(types.AttractionLink)])
  index_mapping = dict(zip(range(len(ppty)), ppty))

  property_df = pd.DataFrame([], columns=["patient_ID"] + list(index_mapping.keys()))
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
              # with default p = 1.5 and T = 2
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
    if normalization == "scaling":
      for column_name in index_mapping.keys():
        property_df[column_name] = property_df[column_name] / max(property_df[column_name])
    elif normalization == "standard":
      for column_name in index_mapping.keys():
        column_mean = property_df[column_name].mean()
        column_std = property_df[column_name].std()
        if column_std == 0:
          property_df = property_df.drop([column_name], axis=1)
        else:
          property_df[column_name] = (property_df[column_name] - column_mean) / column_std
    else:
      property_df = quantile_normalize(property_df)

  # Dump the property vector to CSV (before applying kpca)
  if not normalization:
    normalization = "notnormalized"
  property_vector_pickle_beforekpca = os.path.join(data_dir, 
            "property_vector_beforekpca_p={},T={}_{}_{}.csv".format(p, T, normalization,str(date.today())))
  property_df.to_csv(property_vector_pickle_beforekpca, sep="\t", index=False)

  return property_df

def do_kpca(property_vectors_df):
  def kernel_func(X):
    len_x = X.shape[0]
    dist_array = numpy.random.random((len_x, len_x)).astype(numpy.float16) * 0 - 1.0
    mat = []
    i = 0
    cnt = 0
    total = len_x ** 2
    for a in X:
      a = a.toarray().flatten()
      row = []
      j = 0
      print("--- Working on {}/{} ...".format(cnt, total))
      for b in X:
        b = b.toarray().flatten()
        cnt += 1
        if 0 <= dist_array[i, j]:
          row.append(dist_array[i, j])
        elif 0 <= dist_array[j, i]:
          row.append(dist_array[j, i])
        else:
          # dist = fuzzy_jaccard(a, b)
          dist = tanimoto(a, b)
          row.append(dist)
          dist_array[i, j] = dist
        j += 1
      mat.append(numpy.array(row).astype(numpy.float16))
      i += 1
    return numpy.array(mat, dtype=numpy.float16)

  print("--- Doing KPCA")

  # Compress the vector dataframe to dict
  property_vectors_dict = {}
  for i,j in enumerate(property_vectors_df["patient_ID"]):
    property_vectors_dict[j] = sparse.csr_matrix(property_vectors_df.loc[i,property_vectors_df.columns.drop("patient_ID")].values.tolist())

  X = sparse.vstack(property_vectors_dict.values())
  X_kpca = KernelPCA(kernel = "precomputed").fit_transform(kernel_func(X))

  for k, kpca_v in zip(property_vectors_dict.keys(), X_kpca):
    property_vectors_dict[k] = kpca_v
  return property_vectors_dict

def export_property_vectors(data_dir, property_vectors,p,T, normalization):
  if not normalization:
    normalization = "notnormalized"
  output_file = os.path.join(data_dir, 
          "property_vector_afterkpca_p={},T={}_{}_{}.csv".format(p,T,normalization,str(date.today())))

  print("--- Exporting property vectors to \"{}\"".format(output_file))

  df = dict_to_csv(property_vectors)
  df.to_csv(output_file, sep="\t", index=False)

  if len(zerovector) > 0:
    print(len(zerovector))
    with open(os.path.join(data_dir , "zerovector.txt"), "w") as z:
      z.write("\n".join(zerovector))
  return df

def tanimoto(v1, v2):
  v1_v2 = numpy.dot(v1, v2)
  v1_sq = numpy.sum(numpy.square(v1))
  v2_sq = numpy.sum(numpy.square(v2))
  return v1_v2 / (v1_sq + v2_sq - v1_v2)
