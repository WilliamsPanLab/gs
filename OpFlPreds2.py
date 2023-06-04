import os
import sklearn
import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.model_selection import cross_val_score, train_test_split, GridSearchCV

# Load the g scores
g_df = pd.read_csv('~/g_df_bv.csv')


# Read in the list of subjects that have PropFeats available
prop_feats_dir = '/oak/stanford/groups/leanew1/users/apines/data/gp/PropFeats'
subject_ids = []
for subject_dir in os.listdir(prop_feats_dir):
    if os.path.isdir(os.path.join(prop_feats_dir, subject_dir)):
        subject_ids.append(subject_dir)

# Read in the QC features for each subject and filter out those who do not pass QC
qc_feats_dir = '/oak/stanford/groups/leanew1/users/apines/data/gp/QC_Feats'
remaining_TRs_threshold = 1200

for subject_id in subject_ids:
    qc_feats_file = os.path.join(qc_feats_dir, subject_id, 'QC_Feats.csv')
    qc_feats_df = pd.read_csv(qc_feats_file)
    if qc_feats_df['Var1'][11] < remaining_TRs_threshold:
        subject_ids.remove(subject_id)

# create a dictionary to store the feature data for each subject
subject_features = {}

# loop over the subject ids
for subject_id in subject_ids:
    # define the path to the subject's PropFeats csv
    subject_file = os.path.join(prop_feats_dir, subject_id, 'Prop_Feats.csv')    
    # read in the csv and extract the relevant feature column
    subject_data = pd.read_csv(subject_file, header=None)
    subject_features[subject_id] = subject_data.iloc[:, 1].values.ravel()

# stack the subject vectors on top of each other and calculate column sums
stacked_features = np.vstack(list(subject_features.values()))

# create a DataFrame from the stacked features
df_features = pd.DataFrame(stacked_features, columns=[f'feat_{i}' for i in range(stacked_features.shape[1])])

# insert subject ID as the first column
df_features.insert(0, 'subjectkey', subject_ids)

df = df_features.merge(g_df, on='subjectkey')
# drop subject ids, var1
X=df.drop(df.columns[0],axis=1)
X=X.drop(X.columns[0],axis=1)
# drop var 1
# get y as outcome
Y=df['g']
# drop g from features
X=X.drop(X.columns[-1],axis=1)

# split into train/test sets
X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.3, random_state=42)

# define alphas
alphas = np.exp2(np.arange(16)-10)
# fit ridge cv
ridgeModel = sklearn.linear_model.RidgeCV(alphas=alphas, store_cv_values=True).fit(X_train,y_train)

