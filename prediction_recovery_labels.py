#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 12:58:07 2023

@author: blucie
"""

###############################################################################
# Load libraries
###############################################################################

import pandas as pd
import numpy as np

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import GridSearchCV
from sklearn.pipeline import Pipeline

from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier
from lightgbm import LGBMClassifier

from sklearn.metrics import precision_recall_curve
from sklearn.metrics import balanced_accuracy_score
from sklearn import metrics

import matplotlib.pyplot as plt

###############################################################################
# Functions
###############################################################################

def ohe(data, column_name):
    """ One-hot encoding function for categorical variables
    Parameters
    ----------
    data : pd.DataFrame
        Data frame containing the column to one-hot encode.
        
    column_name : str
        Column name, column should contain a categorical variable.

    Returns
    -------
    data : pd.DataFrame
        Initial data frame, with one-hot encoded column added.
    """
    
    _ohe_cols = pd.get_dummies(data[column_name], dummy_na=False)
    _ohe_cols.columns = [str(col) + '_' + column_name for col in 
                         _ohe_cols.columns]
    data = data.join(_ohe_cols)
    return data

###############################################################################

def calculate_metrics(y_true, y_pred, y_pred_proba, prefix=None):
    """Calculate performance metrics from scores and predictions.

    This is a convenience function for computing performance metrics
    from scores and predictions of the data. Currently the following
    metrics will be calculated:

    - accuracy
    - AUROC
    - AUPRC

    Parameters
    ----------
    y_true : `numpy.array`
        True labels

    y_pred : `numpy.array`
        Predicted labels

    prefix : str, optional
        If set, adds a prefix to the calculated metrics of the form
        'prefix_NAME', where 'NAME' is the name of the metric.

    Returns
    -------
    Dictionary whose keys are the names of the calculated metrics, with
    an optional prefix, and whose values are the respective performance
    measures.
    """
    accuracy = balanced_accuracy_score(y_true, y_pred)
    auprc = metrics.average_precision_score(y_true, y_pred_proba,
                                            pos_label='Recovery')
    auroc = metrics.roc_auc_score(y_true, y_pred_proba)
    fpr, tpr, _ = metrics.roc_curve(y_true, y_pred_proba, 
                                    pos_label='Recovery')
    precision, recall, _ = precision_recall_curve(y_true, y_pred_proba, 
                                                  pos_label='Recovery')

    kv = [
        ('accuracy', accuracy),
        ('auprc', auprc),
        ('auroc', auroc),
    ]

    return dict(kv), fpr, tpr, precision, recall
    
###############################################################################

def preprocessing(time_window):
    """ Preprocessing of the data including one-hot encoding, normalising
    and X/y split
    Parameters
    ----------
    time_window : str 
        Correspond to the weeks to consider to define cohort with very severe
        initial injury.
        Possible values are '1', '4', '1.4'.

    Returns
    -------
    X: pd.DataFrame
        A dataframe with all predictors to use for the prediction task
        
    y: vector
        A vector with the target label (binary)
    """    
    
    # Import data
    df_window_1 = pd.read_csv('/Users/blucie/Desktop/git-repos/outstanding-recovery-SCI/outputs/df_window.' + time_window + '.csv')
    data = df_window_1.copy()

    ## Pre-processing
    
    # Select for variables of interest
    df_subset = data[['SEXCD', 'AGE', 'Lower01', 'Upper01', 'level', 'TX1_R',
              "CORDCD01","CORDOC01", "DECOCD01", "DECOMC01", 
              "DURACD01", "DUROCD01", "EXCDCD01", "EXCICD01", "EXDICD01",
              "FUSCD01", "FUSICD01", "IMMOBC01", "INFXCD01", "INOPCD01", "INTFCD01",
              "INTRAC01", "LAMCD101", "time_decompression_surgery", 'recover', 'new_PTID']]

    # Binarise time to decompression using 24h threshold
    df_subset['time_decompression_surgery'] = df_subset['time_decompression_surgery'].replace([np.nan], [-1])
    df_subset['time_decompression_bin'] = pd.cut(df_subset.time_decompression_surgery, 
                                                 bins=[0, 1440, 15400], labels=[1, 0])
    
    # For categorical data: One-hot-encoding
    cat_cols = ['level', 'TX1_R', "CORDCD01",
                "CORDOC01", "DECOCD01", "DECOMC01", "DURACD01", "DUROCD01", 
                "EXCDCD01", "EXCICD01", "EXDICD01", "FUSCD01", "FUSICD01", 
                "IMMOBC01", "INFXCD01", "INOPCD01", "INTFCD01", "INTRAC01", 
                "LAMCD101", 'time_decompression_bin', 'SEXCD']
    
    for col in cat_cols:
        df_subset[col] = df_subset[col].replace([9], [-1])
        df_subset[col] = df_subset[col].replace([np.nan], [-1])
        #print(df_subset[col].value_counts(dropna=False))
        if df_subset[col].value_counts(dropna=False).shape[0] > 2:
            df_subset = ohe(df_subset, col)
            
    # For continuous data: standardise
    cont_cols = ['AGE', 'Lower01', 'Upper01']
    features = df_subset[cont_cols].copy()
    scaler = StandardScaler().fit(features.values)
    features = scaler.transform(features.values)
    df_subset[cont_cols] = features
    
    #print(df_subset.columns)
    df_subset['SEXCD'] = df_subset['SEXCD'].replace(['1','2'], ['0', '1'])
    df_subset['level'] = df_subset['level'].replace(['C','T'], ['0', '1'])
    
    df_subset_prediction =  df_subset.drop(columns=['time_decompression_surgery',
                                'level', 'TX1_R', "CORDCD01",
                                "CORDOC01", "DECOCD01", "DECOMC01", "DURACD01", "DUROCD01", 
                                "EXCDCD01", "EXCICD01", "EXDICD01", "FUSCD01", "FUSICD01", 
                                "IMMOBC01", "INFXCD01", "INOPCD01", "INTFCD01", "INTRAC01", 
                                "LAMCD101", 'time_decompression_bin', 'SEXCD']).copy()
    df_subset_prediction = df_subset_prediction.dropna()
    
    print(df_subset_prediction.recover.value_counts())
    
    ## Prediction task
    
    # Split between features and labels
    y = df_subset_prediction[['recover']]
    X = df_subset_prediction.drop(columns=['recover', 'new_PTID']).copy()
    
    return X, y

###############################################################################

def preprocessing_baseline(time_window): 
    """ Preprocessing of the data including one-hot encoding, normalising
    and X/y split
    Difference with the previous function: only include baseline info, no
    surgery predictors are included in X when using this function
    Parameters
    ----------
    time_window : str 
        Correspond to the weeks to consider to define cohort with very severe
        initial injury.
        Possible values are '1', '4', '1.4'.

    Returns
    -------
    X: pd.DataFrame
        A dataframe with all predictors to use for the prediction task
        
    y: vector
        A vector with the target label (binary)

    """
    
    # Import data
    df_window_1 = pd.read_csv('/Users/blucie/Desktop/git-repos/outstanding-recovery-SCI/outputs/df_window.' + time_window + '.csv')
    data = df_window_1.copy()

    ## Pre-processing
    
    # Select for variables of interest
    df_subset = data[['SEXCD', 'AGE', 'Lower01', 'Upper01', 'level', 
               'TX1_R', 'recover', 'new_PTID']]
    
    # For categorical data: One-hot-encoding
    cat_cols = ['level', 'TX1_R', 'SEXCD']
    
    for col in cat_cols:
        df_subset[col] = df_subset[col].replace([9], [-1])
        df_subset[col] = df_subset[col].replace([np.nan], [-1])
        #print(df_subset[col].value_counts(dropna=False))
        if df_subset[col].value_counts(dropna=False).shape[0] >= 2:
            df_subset = ohe(df_subset, col)
            
    # For continuous data: standardise
    cont_cols = ['AGE', 'Lower01', 'Upper01']
    features = df_subset[cont_cols].copy()
    scaler = StandardScaler().fit(features.values)
    features = scaler.transform(features.values)
    df_subset[cont_cols] = features
    
    #print(df_subset.columns)
    
    df_subset_prediction =  df_subset.drop(columns=[
        'level', 'TX1_R', 'SEXCD', 'T_level', '2_SEXCD']).copy()
    df_subset_prediction = df_subset_prediction.dropna()
    
    print(df_subset_prediction.recover.value_counts())
    
    ## Prediction task
    
    # Split between features and labels
    y = df_subset_prediction[['recover']]
    X = df_subset_prediction.drop(columns=['recover', 'new_PTID']).copy()
    
    return X, y

###############################################################################
 
def fit_model (X, y, model, param_grid):
    """ Model fitting
    Parameters
    ----------
    X: pd.DataFrame
        A dataframe with all predictors to use for the prediction task
        
    y: vector
        A vector with the target label (binary)
        
    model: 
        Model for fitting and predicting
        
    param_grid: dictionary
        Dictionary including all hyperparameters and potential values for
        tuning in the cross validation step

    Returns
    -------
    df_results: pd.DataFrame
        DataFrame with output from calculate_metrics function
        Includes balanced accuracy, AUROC and AUPRC across 10 random seeds
        
    tprs: list
        List of true positive rate value accross 10 different random seeds 
        Should be of length 10
        
    recalls: list
        List of recall value accross 10 different random seeds 
        Should be of length 10
        
    no_skills:
        Prevalence of the minority class, i.e. performance expected if the
        model did not learn anything from the data
        
    df_feat_importances: pd.DataFrame
        Dataframe with 2 column, one is the feature name and the second is its
        corresponding weight/beta coefficient after fitting the logistic
        regression
        
    """
    
    #Repeat the process with 10 different seeds
    
    df_results = pd.DataFrame({'metrics':['accuracy', 'auprc', 'auroc']})
    tprs = []
    recalls = []
    no_skills = []
    #df_feat_importances = pd.DataFrame()
    base_fpr = np.linspace(0, 1, 101)
    
    for seed in [344, 172, 188, 270, 35, 164, 545, 480, 89, 409]:
        #split the dataset into training (70%) and testing (30%) sets
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.3, random_state = seed, stratify = y)
        
        # Create cross-validation pipeline
        pipeline = Pipeline(steps=[('scaler', None), ('model', model)])

        grid_search = GridSearchCV(
            pipeline,
            param_grid=param_grid,
            cv=5,
            scoring='balanced_accuracy',
            n_jobs=-1,
            verbose=1,
            refit='balanced_accuracy'
        )
        
        # Fit best model
        grid_search.fit(X_train, y_train.values.ravel())

        y_pred_proba = grid_search.predict_proba(X_test)[::,1]
        y_pred = grid_search.predict(X_test)
        
        dict_metrics, fpr, tpr, precision, recall = calculate_metrics(y_test, y_pred, y_pred_proba, prefix='_seed'+str(seed))
        tpr = np.interp(base_fpr, fpr, tpr)
        recall = np.interp(base_fpr, precision, recall)
        
        no_skill = y_test.value_counts()[1] / (y_test.value_counts()[0] + y_test.value_counts()[1])
        
        tprs.append(tpr)
        recalls.append(recall)
        no_skills.append(no_skill)
        
        #importance = log_regression.coef_[0]
        #feat_importances = pd.Series(importance, index = X_test.columns)
        #df_feat_importances = df_feat_importances.append(feat_importances, ignore_index=True)
        
        df_results['seed'+str(seed)] = df_results['metrics'].map(dict_metrics)
    
    return df_results, tprs, recalls, no_skills#, df_feat_importances

###############################################################################
    
def output_metrics (df_results):
    """ Compute meand and SD of the metrics across 10 random seeds
    Parameters
    ----------
    df_results: pd.DataFrame
        DataFrame with output from calculate_metrics function
        Includes balanced accuracy, AUROC and AUPRC across 10 random seeds

    Returns
    -------
    dict(kv): dictionary
        Dictionary containing the average and standard deviation of each metric
        across the 10 random seeds
        
    """
    
    avg_accuracy = df_results.loc[:,[c for c in df_results.columns if c!= "metrics"]].iloc[:1].mean(axis=1)[0]
    avg_auprc = df_results.loc[:,[c for c in df_results.columns if c!= "metrics"]].iloc[1:2].mean(axis=1)[1]
    avg_auroc = df_results.loc[:,[c for c in df_results.columns if c!= "metrics"]].iloc[2:3].mean(axis=1)[2]
    
    std_accuracy = df_results.loc[:,[c for c in df_results.columns if c!= "metrics"]].iloc[:1].std(axis=1)[0]
    std_auprc = df_results.loc[:,[c for c in df_results.columns if c!= "metrics"]].iloc[1:2].std(axis=1)[1]
    std_auroc = df_results.loc[:,[c for c in df_results.columns if c!= "metrics"]].iloc[2:3].std(axis=1)[2]
    
    kv = [
        ('accuracy (mean)', avg_accuracy),
        ('auprc (mean)', avg_auprc),
        ('auroc (mean)', avg_auroc),
        ('accuracy (std)', std_accuracy),
        ('auprc (std)', std_auprc),
        ('auroc (std)', std_auroc),
        
    ]
  
    return dict(kv)

###############################################################################

def roc_curve_plot (tprs):
    """ Plot ROC curve
    Parameters
    ----------
    tprs: list
        List of true positive rate value accross 10 different random seeds 
        Should be of length 10

    Returns
    -------
    fig: matplotlib figure
        Plot of the ROC curve including the expected performance of a random
        classifier
    """
    base_fpr = np.linspace(0, 1, 101)
    
    tprs = np.array(tprs)
    mean_tprs = tprs.mean(axis=0)
    std_tprs = tprs.std(axis=0)
    
    tprs_upper = np.minimum(mean_tprs + std_tprs, 1)
    tprs_lower = mean_tprs - std_tprs
    
    fig = plt.figure(figsize=(5, 5))
    plt.axes().set_aspect('equal', 'datalim')
    plt.plot(base_fpr, mean_tprs, 'b')
    plt.fill_between(base_fpr, tprs_lower, tprs_upper, color='grey', alpha=0.3)
    
    plt.plot([0, 1], [0, 1],'r--')
    plt.xlim([-0.01, 1.01])
    plt.ylim([-0.01, 1.01])
    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive Rate')
    
    return fig

###############################################################################

def pr_curve_plot (recalls, no_skills):
    """ Plot precision-recall curve
    Parameters
    ----------
    recalls: list
        List of recall value accross 10 different random seeds 
        Should be of length 10
        
    no_skills:
        Prevalence of the minority class, i.e. performance expected if the
        model did not learn anything from the data

    Returns
    -------
    fig: matplotlib figure
        Plot of the precision-recall curve including the expected performance 
        of a random classifier
    """
    base_fpr = np.linspace(0, 1, 101)
    recalls = np.array(recalls)
    mean_recalls = recalls.mean(axis=0)
    std_recalls = recalls.std(axis=0)
    
    no_skills = np.array(no_skills)
    mean_no_skills = no_skills.mean(axis=0)
    
    recalls_upper = np.minimum(mean_recalls + std_recalls, 1)
    recalls_lower = mean_recalls - std_recalls
    
    fig = plt.figure(figsize=(5, 5))
    plt.axes().set_aspect('equal', 'datalim')
    plt.plot(base_fpr, mean_recalls, 'b')
    plt.fill_between(base_fpr, recalls_lower, recalls_upper, color='grey', alpha=0.3)
    
    plt.plot([0, 1], [mean_no_skills, mean_no_skills], 'r--')
    plt.xlim([-0.01, 1.01])
    plt.ylim([-0.01, 1.01])
    plt.ylabel("Precision")
    plt.xlabel("Recall")
    
    return fig

###############################################################################

def feat_imp_plot (df_feat_importances, c):
    """ Plot feature importance for the logistic regression
    It will plot the 10 predictors with largest weight when fitting a logistic 
    regression, in decreasing order (less weight on top to most weight at the
    bottom). x-axis is the weight (i.e. beta coefficient), y-axis is the 
    corresponding predictor/feature
    Parameters
    ----------
    df_feat_importances: pd.DataFrame
        Dataframe with 2 column, one is the feature name and the second is its
        corresponding weight/beta coefficient after fitting the logistic
        regression
        
    c: str
        Color to give to the bars

    Returns
    -------
    fig: matplotlib figure
        Plot of the weight of the top 10 features used for the logistic
        regression
    """
    #print(df_feat_importances)
    mean_feat_importance = df_feat_importances.mean()
    idx = mean_feat_importance.abs().nlargest(10).index
    #mean_feat_importance[idx].plot(kind='barh',title = 'Feature Importance')
    fig = plt.figure()
    ax = fig.add_axes([0,0,1,1])
    ax.barh(list(mean_feat_importance[idx].index), list(mean_feat_importance[idx]), color = c)
    
    return fig

###############################################################################


###############################################################################
# Run prediction
###############################################################################

X_1, y_1 = preprocessing('1')
df_results_1, tprs_1, recalls_1, no_skills_1, df_feat_importances_1 = fit_model(X_1, y_1)
output_metrics_1 = output_metrics(df_results_1)
fig_auc_1 = roc_curve_plot(tprs_1)
fig_pr_1 = pr_curve_plot(recalls_1, no_skills_1)
fig_feat_1 = feat_imp_plot(df_feat_importances_1, c='#B3C5E3')

X_4, y_4 = preprocessing('4')
df_results_4, tprs_4, recalls_4, no_skills_4, df_feat_importances_4 = fit_model(X_4, y_4)
output_metrics_4 = output_metrics(df_results_4)
fig_auc_4 = roc_curve_plot(tprs_4)
fig_pr_4 = pr_curve_plot(recalls_4, no_skills_4)
fig_feat_4 = feat_imp_plot(df_feat_importances_4, c='#102E57')

X_1_4, y_1_4 = preprocessing('1.4')
df_results_1_4, tprs_1_4, recalls_1_4, no_skills_1_4, df_feat_importances_1_4 = fit_model(X_1_4, y_1_4)
output_metrics_1_4 = output_metrics(df_results_1_4)
fig_auc_1_4 = roc_curve_plot(tprs_1_4)
fig_pr_1_4 = pr_curve_plot(recalls_1_4, no_skills_1_4)
fig_feat_1_4 = feat_imp_plot(df_feat_importances_1_4, c='#679AE2')

## Baseline - no surgery information
X_1_4_baseline, y_1_4_baseline = preprocessing_baseline('1.4')
df_results_1_4_baseline, tprs_1_4_baseline, recalls_1_4_baseline, no_skills_1_4_baseline, df_feat_importances_1_4_baseline = fit_model(X_1_4_baseline, y_1_4_baseline)
output_metrics_1_4_baseline = output_metrics(df_results_1_4_baseline)
fig_auc_1_4_baseline = roc_curve_plot(tprs_1_4_baseline)
fig_pr_1_4_baseline = pr_curve_plot(recalls_1_4_baseline, no_skills_1_4_baseline)
fig_feat_1_4_baseline = feat_imp_plot(df_feat_importances_1_4_baseline, c='#679AE2')

X_1_baseline, y_1_baseline = preprocessing_baseline('1')

## 

#for model_name in ['lr', 'knn', 'RF']: 
for model_name in ['lightgbm']:
# Define the parameters grid for grid search cross validation
# Define the model used for imputation
    if model_name == 'knn':
        param_grid = [{
            'model__n_neighbors': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
            'model__weights': ['uniform', 'distance']
            }]
        model = KNeighborsClassifier(algorithm='auto', metric='minkowski', p=2)


    elif model_name == 'lr':
        model = LogisticRegression()
        
        param_grid = [{
            'scaler': ['passthrough', StandardScaler()]
            }]

    elif model_name == 'RF':
        param_grid = [{
            'model__bootstrap': [True],
            'model__oob_score': [True, False],
            'model__n_estimators': [15, 25, 50, 75, 100],
            'model__max_features': ['auto', 'sqrt', 'log2'],
            'model__criterion': ['gini', 'entropy'],
            'model__class_weight': [None, 'balanced', 'balanced_subsample']
        }]
        model = RandomForestClassifier(random_state=9)
    
    elif model == 'lightgbm':
        lightgbm = LGBMClassifier(random_state=9)
        
        param_grid = {
            'lightgbm__boosting_type': ['gbdt', 'dart', 'goss', 'rf'],
            'lightgbm__n_estimators': [10, 25, 50, 100, 200],
            'lightgbm__learning_rate': 10.0 ** np.arange(-5, 4),
        }

    df_results_1_4_baseline, tprs_1_4_baseline, recalls_1_4_baseline, no_skills_1_4_baseline = fit_model(X_1_4, y_1_4, model, param_grid)
    output_metrics_1_4_baseline = output_metrics(df_results_1_4_baseline)
    print(model_name)
    print(output_metrics_1_4_baseline)

'''## Baseline week 1-4
# LR
{'accuracy (mean)': 0.5,
 'auprc (mean)': 0.07032406483666331,
 'auroc (mean)': 0.6455284552845529,
 'accuracy (std)': 0.0,
 'auprc (std)': 0.01978154274912521,
 'auroc (std)': 0.08474632161965716}
#knn
{'accuracy (mean)': 0.5, 
 'auprc (mean)': 0.035294117647058816, 
 'auroc (mean)': 0.46524390243902436, 
 'accuracy (std)': 0.0, 
 'auprc (std)': 7.314236392868128e-18, 
 'auroc (std)': 0.016272860895323645}
#RF
{'accuracy (mean)': 0.4981707317073171, 
 'auprc (mean)': 0.0412351161200266, 
 'auroc (mean)': 0.43760162601626007, 
 'accuracy (std)': 0.0029454017776807825, 
 'auprc (std)': 0.011382735093237437, 
 'auroc (std)': 0.11807701122717489}
#lightgbm
{'accuracy (mean)': 0.4981707317073171, 
 'auprc (mean)': 0.0412351161200266, 
 'auroc (mean)': 0.43760162601626007, 
 'accuracy (std)': 0.0029454017776807825, 
 'auprc (std)': 0.011382735093237437, 
 'auroc (std)': 0.11807701122717489}

## week 1-4 (baseline + surgery)
# LR
{'accuracy (mean)': 0.5, 
 'auprc (mean)': 0.159969137108607, 
 'auroc (mean)': 0.6752032520325202, 
 'accuracy (std)': 0.0, 
 'auprc (std)': 0.0920749492997034, 
 'auroc (std)': 0.12219863489881619}
#knn
{'accuracy (mean)': 0.5693089430894309, 
 'auprc (mean)': 0.10080298786181138, 
 'auroc (mean)': 0.608130081300813, 
 'accuracy (std)': 0.0813823552234606, 
 'auprc (std)': 0.06088726178627969, 
 'auroc (std)': 0.10469529535037021}
# RF
{'accuracy (mean)': 0.5463414634146341, 
 'auprc (mean)': 0.17497265286201744, 
 'auroc (mean)': 0.6052845528455284, 
 'accuracy (std)': 0.08188475926022784, 
 'auprc (std)': 0.14919755832512738, 
 'auroc (std)': 0.15890906605445854}
#lightgbm
{'accuracy (mean)': 0.5463414634146341, 
 'auprc (mean)': 0.17497265286201744, 
 'auroc (mean)': 0.6052845528455284, 
 'accuracy (std)': 0.08188475926022784, 
 'auprc (std)': 0.14919755832512738,
 auroc (std)': 0.15890906605445854}


## Baseline week1
#LR
{'accuracy (mean)': 0.5333333333333334, 
 'auprc (mean)': 0.23026586340921015, 
 'auroc (mean)': 0.6365291262135921, 
 'accuracy (std)': 0.03828902431136174, 
 'auprc (std)': 0.09621110435616138, 
 'auroc (std)': 0.08977995625031464}
#knn
{'accuracy (mean)': 0.4983414239482201, 
 'auprc (mean)': 0.11932778893108424, 
 'auroc (mean)': 0.5036003236245955, 
 'accuracy (std)': 0.007146582604554636, 
 'auprc (std)': 0.02818254580333677, 
 'auroc (std)': 0.04573700679406898}
#RF
{'accuracy (mean)': 0.49279935275080905, 
 'auprc (mean)': 0.13791906533852563, 
 'auroc (mean)': 0.5616504854368933, 
 'accuracy (std)': 0.02344014332577166, 
 'auprc (std)': 0.037951520744000725, 
 'auroc (std)': 0.07941249797208955}


## week 1 (baseline + surgery)
# LR
{'accuracy (mean)': 0.5326456310679611, 
 'auprc (mean)': 0.2790324696184552, 
 'auroc (mean)': 0.7302993527508091, 
 'accuracy (std)': 0.04564930247009742, 
 'auprc (std)': 0.10645288707454804, 
 'auroc (std)': 0.05841888660730191}
#knn
{'accuracy (mean)': 0.4993932038834951, 
 'auprc (mean)': 0.14373484553702282, 
 'auroc (mean)': 0.5419093851132687, 
 'accuracy (std)': 0.026499594288728456, 
 'auprc (std)': 0.04876089580254087, 
 'auroc (std)': 0.06863212437929878}
#RF
{'accuracy (mean)': 0.5010517799352752, 
 'auprc (mean)': 0.16821691945953324, 
 'auroc (mean)': 0.5942152103559871, 
 'accuracy (std)': 0.019251205899742652, 
 'auprc (std)': 0.055945948411547035, 
 'auroc (std)': 0.05906655278688844}'''



