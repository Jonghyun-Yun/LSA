import numpy as np
from matplotlib import pyplot as plt

from sklearn.cluster import SpectralCoclustering
from sklearn.metrics import consensus_score

import pandas as pd
D = pd.read_csv('chessB_pn/rbf_dist_zw.csv', sep=',',header=0)
data = D.values

plt.matshow(data, cmap=plt.cm.Blues)
plt.title("Original dataset")

# row_idx = np.linspace(0,(data.shape[0]-1), num=data.shape[0]).astype(int)
# col_idx = np.linspace(0,(data.shape[1]-1), num=data.shape[1]).astype(int)
K = 2

model = SpectralCoclustering(n_clusters=K, random_state=0)
model.fit(data)

fit_data = data[np.argsort(model.row_labels_)]
fit_data = fit_data[:, np.argsort(model.column_labels_)]

plt.matshow(fit_data, cmap=plt.cm.Blues)
plt.title("After biclustering; rearranged to show biclusters")

df_row = pd.DataFrame(model.row_labels_)
df_col = pd.DataFrame(model.column_labels_)

# save the dataframe as a csv file
df_row.to_csv("chessB_pn/cl_row" + str(K) + ".csv")
df_col.to_csv("chessB_pn/cl_col" + str(K) + ".csv")
