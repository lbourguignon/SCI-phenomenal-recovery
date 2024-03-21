#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 15:02:10 2023

@author: blucie
"""

base_fpr = np.linspace(0, 1, 101)
fig = plt.figure(figsize=(5, 5))

tprs_1 = np.array(tprs_1)
mean_tprs_1 = tprs_1.mean(axis=0)

tprs_4 = np.array(tprs_4)
mean_tprs_4 = tprs_4.mean(axis=0)

tprs_1_4 = np.array(tprs_1_4)
mean_tprs_1_4 = tprs_1_4.mean(axis=0)

plt.axes().set_aspect('equal', 'datalim')
plt.plot(base_fpr, mean_tprs_1, '#B3C5E3', label='Week 1')
plt.plot(base_fpr, mean_tprs_4, '#102E57', label='Week 4')
plt.plot(base_fpr, mean_tprs_1_4, '#679AE2', label='Weeks 1 and 4')

plt.plot([0, 1], [0, 1],'r--', label='Random classifier')
plt.xlim([-0.01, 1.01])
plt.ylim([-0.01, 1.01])
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')
plt.legend(loc='lower right')
plt.show()



base_fpr = np.linspace(0, 1, 101)
recalls_1 = np.array(recalls_1)
mean_recalls_1 = recalls_1.mean(axis=0)
no_skills_1 = np.array(no_skills_1)
mean_no_skills_1 = no_skills_1.mean(axis=0)

recalls_4 = np.array(recalls_4)
mean_recalls_4 = recalls_4.mean(axis=0)
no_skills_4 = np.array(no_skills_4)
mean_no_skills_4 = no_skills_4.mean(axis=0)

recalls_1_4 = np.array(recalls_1_4)
mean_recalls_1_4 = recalls_1_4.mean(axis=0)
no_skills_1_4 = np.array(no_skills_1_4)
mean_no_skills_1_4 = no_skills_1_4.mean(axis=0)


fig = plt.figure(figsize=(5, 5))
plt.axes().set_aspect('equal', 'datalim')
plt.plot(base_fpr, mean_recalls_1, '#B3C5E3', label='Week 1')
plt.plot(base_fpr, mean_recalls_4, '#102E57', label='Week 4')
plt.plot(base_fpr, mean_recalls_1_4, '#679AE2', label='Weeks 1 and 4')

plt.plot([0, 1], [mean_no_skills_1, mean_no_skills_1], linestyle = '--', 
         color = '#B3C5E3')
plt.plot([0, 1], [mean_no_skills_4, mean_no_skills_4], linestyle = '--', 
         color = '#102E57')
plt.plot([0, 1], [mean_no_skills_1_4, mean_no_skills_1_4], linestyle = '--', 
         color = '#679AE2', label='Corresponding random classifier')
plt.xlim([-0.01, 1.01])
plt.ylim([-0.01, 1.01])
plt.legend(loc='upper right')
plt.ylabel("Precision")
plt.xlabel("Recall")


# fpr, tpr, _ = metrics.roc_curve(y_test,  y_pred_proba, pos_label='Recovery')
# auc = metrics.roc_auc_score(y_test, y_pred_proba)

# #create ROC curve
# plt.plot(fpr,tpr)
# plt.ylabel('True Positive Rate')
# plt.xlabel('False Positive Rate')
# plt.show()

# precision, recall, thresholds = precision_recall_curve(y_test, y_pred_proba, pos_label='Recovery')
# no_skill = y_test.value_counts()[1] / (y_test.value_counts()[0] + y_test.value_counts()[1])
# plt.plot([0, 1], [no_skill, no_skill], linestyle='--', label='No Skill')
# plt.plot(recall, precision)
# plt.ylabel("Precision")
# plt.xlabel("Recall")
# plt.title("Train Precision-Recall curve")
# plt.show()

# importance = log_regression.coef_[0]
# # summarize feature importance
# for i,v in enumerate(importance):
#  print('Feature: %0d, Score: %.5f' % (i,v))

# #importance is a list so you can plot it. 
# feat_importances = pd.Series(importance, index = X_test.columns)
# idx = feat_importances.abs().nlargest(20).index
# feat_importances[idx].plot(kind='barh',title = 'Feature Importance')