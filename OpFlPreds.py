import os
import sklearn
import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.model_selection import cross_val_score, train_test_split, GridSearchCV

# Load the g scores
g_df = pd.read_csv('~/g_df_bv.csv')

# read in brain features for all subjects
data_dir = '/scratch/users/apines/gp/PropFeats'
subject_features = {}
for filename in os.listdir(data_dir):
    if filename.endswith('.csv'):
        subject_id = filename.split('_')[0]
        subject_data = pd.read_csv(os.path.join(data_dir, filename), header=None)
        subject_features[subject_id] = subject_data.values.ravel()

# stack subject vectors on top of each other and calculate column sums
stacked_features = np.vstack(list(subject_features.values()))

# create a DataFrame from the stacked features
df_features = pd.DataFrame(stacked_features, columns=[f'feat_{i}' for i in range(stacked_features.shape[1])])

# insert subject ID as the first column
subject_ids = list(subject_features.keys())
df_features.insert(0, 'subjectkey', subject_ids)

# drop where any values are 0
df_features = df_features.loc[:, ~(df_features == 0).any()]

# merge with g
df = df_features.merge(g_df, on='subjectkey')

# extract x (all brain features, dropping subjectkey and g) and y (g)
X=df.drop(df.columns[0],axis=1)
Y=df['g']
X=X.drop(X.columns[-1],axis=1)

# split into train/test sets
X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.3, random_state=42)

# define alphas
alphas = np.exp2(np.arange(16)+70)
# fit ridge cv
ridgeModel = sklearn.linear_model.RidgeCV(alphas=alphas, store_cv_values=True).fit(X_train,y_train)

