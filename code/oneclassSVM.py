#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 14:14:19 2024

@author: blucie
"""

import pandas as pd
import numpy as np
from sklearn.svm import OneClassSVM
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import StandardScaler
from scipy.sparse import linalg
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

sygen = pd.read_csv('/Volumes/blucie/PhD/1_SCI/10_miraculous-recovery/modified_sygen.csv')

sygen_sub = sygen[['ptid', 'sexcd', 'age', 'lower04', 'upper04', 'ais4', 'lower_locf', 'upper_locf', 'ais52', 'NLI_calc_bin']]

sygen_sub.dtypes

#sygen_sub['sexcd'] = sygen_sub['sexcd'].astype(object)
sygen_sub['ais52'] = sygen_sub['ais52'].astype(object)
sygen_sub['ais4'] = sygen_sub['ais4'].astype(object)
sygen_sub.dtypes

sygen_sub = pd.concat([sygen_sub,
                       pd.get_dummies(sygen_sub['ais4']).astype(int),
                       pd.get_dummies(sygen_sub['NLI_calc_bin']).astype(int)], 
                      axis=1)
sygen_sub = sygen_sub.rename(columns = {'AIS A': 'A4', 'AIS B': 'B4', 'AIS C': 'C4', 'AIS D': 'D4', 'AIS E': 'E4'})
sygen_sub = pd.concat([sygen_sub,
                       pd.get_dummies(sygen_sub['ais52']).astype(int)],
                      axis=1)
X = sygen_sub[['ptid', 'sexcd', 'age', 'lower04', 'upper04', 'lower_locf',
               'upper_locf', 'c', 't','AIS A', 'AIS B', 'AIS C', 'AIS D', 'AIS E', 'A4', 'B4', 'C4', 'D4', 'E4']] #'AIS A', 'AIS B', 'AIS C', 'AIS D', 'AIS E',
X.dropna(axis = 0, how = 'any', inplace = True)

# rescale the data
x_scaled = StandardScaler().fit_transform(X.drop(['ptid'], axis=1))

# reduce the data to 2 dimensions using t-SNE
x_reduced_TSNE = TSNE(n_components=2, random_state=0).fit_transform(x_scaled)

# reduce the data to 2 dimensions using t-SNE
PCA_model = PCA(n_components=5, random_state=0)
x_reduced_PCA = PCA_model.fit_transform(x_scaled)

sum(PCA_model.explained_variance_ratio_)

x_reduced = x_reduced_TSNE

# fit the model to the reduced data
svm = OneClassSVM(kernel='rbf', gamma='scale', nu=0.05)
svm.fit(x_reduced)

# extract the model predictions
x_predicted = svm.predict(x_scaled)


plt.scatter(x_reduced[:,0], x_reduced[:,1], c = svm.predict(x_reduced))
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.show()


# define the meshgrid
x_min, x_max = x_reduced[:, 0].min() - 5, x_reduced[:, 0].max() + 5
y_min, y_max = x_reduced[:, 1].min() - 5, x_reduced[:, 1].max() + 5

x_ = np.linspace(x_min, x_max, 500)
y_ = np.linspace(y_min, y_max, 500)

xx, yy = np.meshgrid(x_, y_)

# evaluate the decision function on the meshgrid
z = svm.decision_function(np.c_[xx.ravel(), yy.ravel()])
z = z.reshape(xx.shape)

# plot the decision function and the reduced data
plt.contourf(xx, yy, z, cmap=plt.cm.PuBu)
a = plt.contour(xx, yy, z, levels=[0], linewidths=2, colors='darkred')
b = plt.scatter(x_reduced[x_predicted == 1, 0], x_reduced[x_predicted == 1, 1], c='white', edgecolors='k')
c = plt.scatter(x_reduced[x_predicted == -1, 0], x_reduced[x_predicted == -1, 1], c='gold', edgecolors='k')
plt.legend([a.collections[0], b, c], ['learned frontier', 'regular observations', 'abnormal observations'], bbox_to_anchor=(1.05, 1))
plt.axis('tight')
plt.show()

unique, counts = np.unique(svm.predict(x_reduced), return_counts=True)
print(np.asarray((unique, counts)).T)

#clf = OneClassSVM(gamma='auto').fit(X)
#clf.predict(X)
#clf.decision_function(X)
#unique, counts = np.unique(svm.predict(x_reduced), return_counts=True)
#print(np.asarray((unique, counts)).T)
#clf.score_samples(X)

d = {'ptid': X['ptid'].to_numpy(), 'svm_outlier': x_predicted}
df = pd.DataFrame(data=d)

df.to_csv('/Volumes/blucie/PhD/1_SCI/10_miraculous-recovery/outliers_week4_full.csv', index=False)
