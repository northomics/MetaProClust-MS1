#!/usr/bin/env python3
import pandas as pd
import numpy as np
from sklearn.impute import KNNImputer
import argparse


usage = "usage: %prog [options] -i filtered_intensity.csv -o knn_imputed_intensity.csv -d diamond.out"

parser = argparse.ArgumentParser(description="quick knn imputation for peptide intensity")
parser.add_argument("-i", "--input", dest="input", required=True, metavar="filtered_intensity.csv", help="input filtered intensity files")
parser.add_argument("-o", "--output", dest="output", required=True, metavar="knn_imputed_intensity.csv", help="output imputed intensity files")
parser.add_argument("-n", "--n_neighbours", dest="n", required=False, default=5, metavar="n_neighbours", help="number of neighbours in KNN impumentation")
args = parser.parse_args()



#feature_file = "../data/ms1_filtered.csv"
feature_file = args.input
output_file = args.output

## import the filtered and library transformed ms1 features
features = pd.read_csv(feature_file, header = 0, index_col = 0)

#print(features.head())
n_neighbours = args.n
imputer = KNNImputer(n_neighbors=n_neighbours, weights="uniform")
features = imputer.fit_transform(features)
#log_features = np.log2(features)

#np.savetxt("../data/ms1_filtered_knn.csv", features, delimiter=",")
np.savetxt(output_file, features, delimiter=",")
