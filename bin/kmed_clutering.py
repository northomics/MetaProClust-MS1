#!/usr/bin/env python3
import pandas as pd
import numpy as np
import os
#import matplotlib.pyplot as plt
from scipy import signal
from sklearn.decomposition import PCA, FastICA
from sklearn.metrics import pairwise_distances, silhouette_score
from sklearn.cluster import AgglomerativeClustering
from sklearn_extra.cluster import KMedoids
import argparse
#from os.path import splitext, basename, join, exists
# scipy.spatial.distance is where the correlation comes from?


#################
### FUNCTIONS ###
#################
# to find the best k number of clusters
#https://github.com/LaureBerti/Learn2Clean/blob/master/python-package/learn2clean/clustering/clusterer.py
def compare_k_AggClustering(k_list, X):
    # Run clustering with different k and check the metrics
    silhouette_list = []
    for p in k_list:
        print("Calculating silhouette score for k =", p)
        clusters = AgglomerativeClustering(n_clusters=p, affinity='precomputed', linkage='average', distance_threshold=None)
        clusters.fit(X)
        # The higher (up to 1) the better
        s = round(silhouette_score(X, clusters.labels_, metric="precomputed"), 4)
        silhouette_list.append(s)
    # The higher (up to 1) the better
    key = silhouette_list.index(max(silhouette_list))
    k = k_list.__getitem__(key)
    print("AggClust best silhouette =", max(silhouette_list), " for k =", k)
    return k, silhouette_list



def compare_k_med(k_list, X):
    # Run clustering with different k and check the metrics
    silhouette_list = []
    inertia_list = []
    for p in k_list:
        print("Calculating silhouette score for k =", p)
        clusters = KMedoids(n_clusters=p, metric='precomputed', random_state=2248, init='k-medoids++')
        clusters.fit(X)
        # The higher (up to 1) the better
        s = round(silhouette_score(X, clusters.labels_, metric="precomputed"), 4)
        silhouette_list.append(s)
        i = clusters.inertia_
        inertia_list.append(i)
    # The higher (up to 1) the better
    key = silhouette_list.index(max(silhouette_list))
    k = k_list.__getitem__(key)
    print("Kmed best silhouette =", max(silhouette_list), " for k =", k)
    return k, silhouette_list, inertia_list

####################################################
### SET UP CMND LINE OPTIONS FOR FILE IN/OUTPUT ###
####################################################

usage = "usage: %prog [options] -s s_mat.csv -n nthreads -d directory/"

parser = argparse.ArgumentParser(description="Peptide/feature clusters from pre-computed ICA decomposition")
parser.add_argument("-s", "--s_mat", dest="s_mat", required=True, metavar="S.csv", help="S matrix ICA output")
parser.add_argument("-d", "--dir", dest="dir", required=True, metavar="directory/", help="directory for your output files")
parser.add_argument("-n", "--threads", dest="n_threads", default=1, metavar="nthreads/", help="number of threads to use")
args = parser.parse_args()



## data input files
s_mat = args.s_mat
output_dir = os.path.dirname(args.dir)

## data output files
kmed_labels = os.path.join(output_dir + "/k-med_labels.csv")
sil_kmed = os.path.join(output_dir + "/k-med_silhouette.csv")
inertia_kmed = os.path.join(output_dir + "/k-med_inertia.csv")


## import the normalized and transformed ms1 features


### ICA
# https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.FastICA.html
# https://stackoverflow.com/questions/35407669/independent-component-analysis-ica-in-python
# using the number of components that explained 99% of the variance using PCA
## following Sastry2019a
#ica = FastICA(random_state=2248, n_components=213, whiten=True, algorithm='parallel',  tol=10e-8)

# Compute ICA...already computed//
#S_ = ica.fit_transform(log_features)  # Reconstruct signals
#A_ = ica.mixing_  # Get estimated mixing matrix

S_ = pd.read_csv(s_mat, header=0, index_col=0)
#A_ = pd.read_csv(a_mat, header=0)

nthreads = int(args.n_threads)

## Compute distance matrix for mixing 
mixing_cor = pairwise_distances(S_, Y=None, metric='correlation', n_jobs=nthreads) ##  increase the njobs if run on beter server

## Trying Kmedoids
## see which k has best silhouette score
k_list = list(range(2, 50))

best_k_med, kmed_sil_list, kmed_inert  = compare_k_med(k_list, mixing_cor)
kmedoids_S = KMedoids(n_clusters=best_k_med, metric='precomputed', random_state=2248, init='k-medoids++').fit(mixing_cor)
kmed_clust_labels = kmedoids_S.labels_

## save these asap so we can start working on this asap
np.savetxt(kmed_labels, kmed_clust_labels.astype(int), delimiter=',', fmt='%i')
np.savetxt(sil_kmed, kmed_sil_list, delimiter=',', fmt='%f')
np.savetxt(inertia_kmed, kmed_inert, delimiter=',', fmt='%f')

