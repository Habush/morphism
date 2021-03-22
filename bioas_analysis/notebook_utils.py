import warnings
import numpy as np
import pandas as pd
warnings.filterwarnings('ignore')
from datetime import datetime
from sklearn.model_selection import RandomizedSearchCV
from sklearn.metrics import roc_auc_score, recall_score, average_precision_score, balanced_accuracy_score, precision_score, make_scorer, confusion_matrix
from sklearn.model_selection import cross_val_score, cross_validate, StratifiedKFold, train_test_split
from xgboost import XGBClassifier
import matplotlib.pyplot as plt
from sklearn.preprocessing import LabelEncoder
from sklearn.decomposition import KernelPCA, PCA

def avp_0(y_true, y_pred):
    avp = average_precision_score(y_true, y_pred, pos_label=0)
    return avp

average_precision_0 = make_scorer(avp_0, greater_is_better=True)

seed = 42
st_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
score_cols = ["test_balanced_accuracy", "test_recall_0","test_precision_0",
               "test_recall_1","test_precision_1", "test_auc"]

def calc_scores(clf, X_test, y_test):
    y_pred = clf.predict(X_test)
    recall_0, recall_1 = recall_score(y_test, y_pred, pos_label=0), recall_score(y_test, y_pred, pos_label=1)
    precision_0, precision_1 =  precision_score(y_test, y_pred, pos_label=0), precision_score(y_test, y_pred, pos_label=1)
    acc = balanced_accuracy_score(y_test, y_pred)
    auc_score = roc_auc_score(y_test, clf.predict_proba(X_test)[:,1])
    tn, fp, fn, tp = confusion_matrix(y_test, y_pred, labels=[0,1]).ravel()
    print("tn: {}, fp:{}, fn:{}, tp:{}" .format(tn, fp, fn, tp))
    print(confusion_matrix(y_test, y_pred, labels=[0,1]))
    fpr = fp/(fp+tn)
    fnr = fn/(tp+fn)
    arr = np.array([[acc, recall_0, precision_0,recall_1, precision_1,auc_score, fpr, fnr]])
    return pd.DataFrame(data=arr, columns=["balanced_accuracy", "recall_0", "precision_0", "recall_1", "precision_1", "auc", "False_positive_rate","False_negative_rate"])

def recall_0(y_true, y_pred):
    return recall_score(y_true, y_pred, pos_label=0)

def precision_0(y_true, y_pred):
    return precision_score(y_true, y_pred, pos_label=0)


def timer(start_time=None):
    if not start_time:
        start_time = datetime.now()
        return start_time

    elif start_time:
        thour, temp_sec = divmod((datetime.now() - start_time).total_seconds(), 3600)
        tmin, tsec = divmod(temp_sec, 60)
        print('\n Time taken: %i hours %i minutes and %s seconds.' % (thour, tmin, round(tsec, 2)))

def param_tuning(X, y, n_folds=5, param_comb=25, scoring='roc_auc', jobs=12, scale_pos_weight=False):
    if scale_pos_weight:
        xgb = XGBClassifier(learning_rate=0.02, n_estimators=600, objective='binary:logistic',
                    silent=True, nthread=1, scale_pos_weight=scale_pos_weight)
    else:
        xgb = XGBClassifier(learning_rate=0.02, n_estimators=600, objective='binary:logistic',
                    silent=True, nthread=1)
    skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=42)
    params = {#'n_estimators': [30, 50, 80, 100, 120],
              'n_estimators': [300, 400, 500, 600, 700],
              'learning_rate': [0.01, 0.02, 0.03, 0.05, 0.07],
              'gamma': [0.5, 1, 1.5, 2, 5],
#               'gamma': [8, 9, 10, 11, 12],
              'max_depth': [3, 4, 5, 6],
              'subsample': [0.6, 0.8, 1.0],
              'colsample_bytree': [0.6, 0.8, 1.0],
              'min_child_weight': [1, 2, 3, 4, 5],
              'max_delta_step': [1,2,3,4,5]
              }
    if scale_pos_weight:
        params['scale_pos_weight'] = [0.1, 0.2, 0.36, 1]
    rand_search = RandomizedSearchCV(xgb, param_distributions=params, n_iter=param_comb, scoring=scoring, n_jobs=jobs,
                                   cv=skf.split(X, y), verbose=3, random_state=42)

    start_time = timer(None) # timing starts from this point for "start_time" variable
    rand_search.fit(X, y)
    timer(start_time)
    print("Best Score: {:.3%}".format(rand_search.best_score_))
    print(rand_search.best_params_)
    return rand_search


def get_scores(cv_results, score_keys=None, df_cols=None):
    if score_keys is None:
        score_keys = score_cols
    if df_cols is None:
        df_cols = score_cols
    scores = np.empty([1, len(score_keys)])
    for i, s in enumerate(score_keys):
        scores[0][i] = np.mean(cv_results[s])
    scores_df = pd.DataFrame(data=scores, columns=[i.replace("test_","") for i in df_cols])
    return scores_df

def generate_X_y(df, dim, train, test):
    X,y = df[df.columns.drop("posOutcome")], df["posOutcome"]
    if dim:
        X = X[X.columns[:dim]]
    if train:
        X_train, X_test, y_train, y_test = X.loc[train], X.loc[test], y.loc[train], y.loc[test]
    else:    
        X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y, test_size=0.3, random_state=42)
    return X_train, X_test, y_train, y_test

def do_train_val_test(df, dim=False, scoring=False, scale_pos_weight=False, train=False, test=False):
    X_train, X_test, y_train, y_test = generate_X_y(df, dim, train, test)
    print("---- Do parameter tuning")
    if scoring:
        print("Scoring parameter {}".format(scoring))
        rand_search = param_tuning(X_train, y_train, scoring=scoring, scale_pos_weight=scale_pos_weight)
    else:
        print("default Scoring parameter will be used")
        rand_search = param_tuning(X_train, y_train)
    clasifier = rand_search.best_estimator_
    scoring_val = {"balanced_accuracy": make_scorer(balanced_accuracy_score),
           "recall_0": make_scorer(recall_0), "precision_0": make_scorer(precision_0),
           "recall_1": make_scorer(recall_score), "precision_1": make_scorer(precision_score), "auc": "roc_auc" }
    
    validation = cross_validate(clasifier, X_train, y_train, scoring=scoring_val,
                                 n_jobs=-1, verbose=1, return_train_score=True, cv=st_cv)
    val_score = get_scores(validation, score_cols)
    print("--- Validation score")
    print(val_score.mean())
    
    # on test dataset

    clasifier.fit(X_train, y_train)
    test_scores = calc_scores(clasifier, X_test, y_test)

    print("--- Test score")
    print(test_scores.mean())
    
    return val_score, test_scores

def print_comparision(df_dict):
    df_dict = dict(sorted(df_dict.items()))
    cols = df_dict.keys()
    rows = []
    for df in df_dict.values():
        rows = rows + df.columns.tolist()
    rows = sorted(set(rows))
    df = pd.DataFrame([], columns=cols, index=set(rows))
    for r in rows:
        df.loc[r] = [d[r].values[0] for d in df_dict.values()]
    df = df.sort_index()
    return df

def find_best_dim(emb_vector, scoring=False, scale_pos_weight=False, train_set=False, test_set=False):
    results = {}
    best_auc = 0
    low_fpr = 1
    best_balanced_acc = 0
    for i in range(2,150):
        print("------------ Dimention {}".format(i))
        val, test = do_train_val_test(emb_vector, dim=i, scoring=scoring, 
                                      scale_pos_weight=scale_pos_weight, train=train_set, test=test_set)
        results["dim {}".format(i)] = {"val":val, "test":test}
        auc = val["auc"].values[0]
        fpr = test["False_positive_rate"].values[0]
        balanced_acc = val["balanced_accuracy"].values[0]
        if  auc > best_auc:
            best_auc = auc
            dim_best_auc = i
        if fpr < low_fpr:
            low_fpr = fpr
            dim_low_fpr = i
        if balanced_acc > best_balanced_acc:
            best_balanced_acc = balanced_acc
            best_acc_dim = i
    results["dim all"] = {"val":val, "test":test}
    print("Best auc= {}, dimentions used= {}".format(best_auc, dim_best_auc))
    print("Best Balanced_accuracy= {}, dimentions used= {}".format(best_balanced_acc, best_acc_dim))
    print("Lowest FPR = {}, dimentions used= {}".format(low_fpr, dim_low_fpr))
    val, test = do_train_val_test(emb_vector, scoring=scoring, scale_pos_weight=scale_pos_weight, train=train_set, test=test_set)
    return results

def do_kpca_rbf(df, test_data=False, vec_space=False, kernel=False):
    if kernel:
        kpca = KernelPCA(kernel=kernel, n_components=150)
    else:
        kpca = KernelPCA(kernel="rbf")
    if test_data:
        kpca.fit(vec_space[vec_space.columns.drop(["patient_ID"])])
        X_kpca = kpca.transform(df[df.columns.drop(["patient_ID"])])
    else:
        X_kpca = kpca.fit_transform(df[df.columns.drop(["patient_ID"])])
    df_res = pd.DataFrame(X_kpca, columns=list(range(X_kpca.shape[1])))
    df_res["patient_ID"] = df["patient_ID"]
    return df_res

def do_pca(df, test_data=False, vec_space=False):
    pca = PCA(n_components=2)
    if test_data:
        pca.fit(vec_space[vec_space.columns.drop(["patient_ID"])])
        X_pca = pca.transform(df[df.columns.drop(["patient_ID"])])
    else:
        X_pca = pca.fit_transform(df[df.columns.drop(["patient_ID"])])
    df_res = pd.DataFrame(X_pca, columns=list(range(X_pca.shape[1])))
    df_res["patient_ID"] = df["patient_ID"]
    return df_res
    