Method: <sup id="b5b03d98105e15b5544548c3e73abbea"><a href="#dhillon_co-clustering_2001" title="Dhillon, Co-Clustering Documents and Words Using Bipartite Spectral Graph Partitioning, 269--274, in in: {Proceedings of the Seventh {{ACM SIGKDD}} International Conference on {{Knowledge}} Discovery and Data Mining}, edited by {Association for Computing Machinery} (2001)">dhillon_co-clustering_2001</a></sup> Dhillon, I. S., Co-clustering documents and words using bipartite spectral graph partitioning, In , Proceedings of the Seventh {{ACM SIGKDD}} International Conference on {{Knowledge}} Discovery and Data Mining (pp. 269–274) (2001). {New York, NY, USA}: {Association for Computing Machinery}.

Code: <https://scikit-learn.org/stable/modules/biclustering.html> <https://scikit-learn.org/stable/auto_examples/bicluster/plot_bicluster_newsgroups.html#sphx-glr-auto-examples-bicluster-plot-bicluster-newsgroups-py>


# initialization

`w0` and `z0` are matrices for latent item and respondent coordinates.

```R
out_dir <- "chessB_pn/"
num_chain <- 3
double_z <- 0
double_w <- 0
HAS_REF <- 0
library(art)
library(coda)
library(dplyr)
library(stringr)
library(magrittr)
library(bayesplot)
library(foreach)
library(doParallel)
## registerDoParallel(cores = detectCores() - 1)
stopImplicitCluster()
registerDoParallel(2)

setwd("/Users/yunj/Dropbox/research/lsjm-art/lsjm-code")

w0 <- readr::read_csv(paste0(out_dir, "w0.csv"))
z0 <- readr::read_csv(paste0(out_dir, "z0.csv"))
```


# calculate Euclidean distance and Gaussian distance over all item-respondent pairs

-   `dist_zw.csv` Euclidean dist
-   `rbf_dist_zw.csv` Gaussian dist

```R
D <- matrix(0, nrow = nrow(z0), ncol = nrow(w0))
for (i in 1:nrow(D)) {
  for (j in 1:ncol(D)) {
    D[i, j] <- sqrt(sum((z0[i, ] - w0[j, ])^2))
  }
}

RBF_D <- exp(-1 / 2 * D)
readr::write_csv(as.data.frame(D), paste0(out_dir, "dist_zw.csv"))
readr::write_csv(as.data.frame(RBF_D), paste0(out_dir, "rbf_dist_zw.csv"))
```


# biclustering using Python

See <https://scikit-learn.org/stable/modules/biclustering.html> for the spectral co-clustering method we implement below: ( <https://scikit-learn.org/stable/auto_examples/bicluster/plot_spectral_coclustering.html#sphx-glr-auto-examples-bicluster-plot-spectral-coclustering-py>)

```jupyter-python
import numpy as np
from matplotlib import pyplot as plt

from sklearn.cluster import SpectralCoclustering
from sklearn.metrics import consensus_score

import pandas as pd
D = pd.read_csv('chessB_pn/rbf_dist_zw.csv', sep=',',header=0)
data = D.values

plt.matshow(data, cmap=plt.cm.Blues)
plt.title("Original dataset")
```

```jupyter-python
# row_idx = np.linspace(0,(data.shape[0]-1), num=data.shape[0]).astype(int)
# col_idx = np.linspace(0,(data.shape[1]-1), num=data.shape[1]).astype(int)
K = 2

model = SpectralCoclustering(n_clusters=K, random_state=0)
model.fit(data)

fit_data = data[np.argsort(model.row_labels_)]
fit_data = fit_data[:, np.argsort(model.column_labels_)]

plt.matshow(fit_data, cmap=plt.cm.Blues)
plt.title("After biclustering; rearranged to show biclusters")
```

```jupyter-python
df_row = pd.DataFrame(model.row_labels_)
df_col = pd.DataFrame(model.column_labels_)

# save the dataframe as a csv file
df_row.to_csv("chessB_pn/cl_row" + str(K) + ".csv")
df_col.to_csv("chessB_pn/cl_col" + str(K) + ".csv")
```


# back to R for visualization

```R
source("R/art-functions.R")

K <- 2

w0 <- readr::read_csv(paste0(out_dir, "w0.csv"))
z0 <- readr::read_csv(paste0(out_dir, "z0.csv"))

cl_z <- readr::read_csv(paste0(out_dir, "cl_row", K, ".csv"))[, 2]
cl_w <- readr::read_csv(paste0(out_dir, "cl_col", K, ".csv"))[, 2]

xmin <- min(z0[, 1], w0[, 1])
ymin <- min(z0[, 2], w0[, 2])
xmax <- max(z0[, 1], w0[, 1])
ymax <- max(z0[, 2], w0[, 2])

pdf(paste0(out_dir, "spc_cl_dist_rbf", K, ".pdf"))

print(cl_lsjmplot(z0, w0, cl_z, cl_w, xlim = c(xmin, xmax), ylim = c(ymin, ymax)))

dev.off(which = dev.cur())
```
