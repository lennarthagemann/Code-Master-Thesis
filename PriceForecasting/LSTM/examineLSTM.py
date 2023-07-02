import os
from tensorflow.keras.models import Sequential, load_model
from tensorflow.keras.layers import Dense, LSTM, Dropout, InputLayer
from tensorflow.keras.callbacks import ModelCheckpoint
from tensorflow.keras.losses import MeanSquaredError
from tensorflow.keras.metrics import RootMeanSquaredError
from tensorflow.keras.optimizers import Adam
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import numpy as np

df_path = r"C:\Users\lenna\OneDrive - NTNU\Master Thesis\Data\Spot Price Data"
files = os.listdir(df_path)

csv_files = [f for f in files if f.endswith(".csv")]
df_list = []
i = 0

for f in csv_files:
    i += 1
    df = pd.read_csv(os.path.join(df_path, f), header=0,
                     sep=";", index_col=0, decimal=',')
    mask = df.columns.str.endswith('volume')  
    df.columns = df.columns.str.replace('Nord Pool ', '')  
    df_list.append(df.loc[:,~mask])

stacked_dfs = pd.DataFrame()
test = pd.Series(df_list[0].T.iloc[0])
test = test.set_axis( pd.to_datetime(test.name) + pd.to_timedelta(test.index - 1, unit='h'))
test.name = 'price'

dfs_to_stack = []
for i in range(len(df_list)):
    for row in df_list[i].T.index:
        df = pd.Series(df_list[i].T.loc[row])
        df = df.set_axis( pd.to_datetime(df.name) + pd.to_timedelta(df.index - 1, unit='h'))
        df.name = 'price'
        dfs_to_stack.append(df)
stacked_dfs = pd.concat(dfs_to_stack)
stacked_dfs.sort_index(inplace=True)
# Get the date of the first index in stacked_dfs, and add an hour for the number of each row
first_date = stacked_dfs.index[0].date()
stacked_dfs = stacked_dfs.set_axis(pd.date_range(start = first_date, periods=len(stacked_dfs), freq='H'))

# Replace nan rows with mean value of surrounding rows
stacked_dfs_nan =  stacked_dfs.isna()
def replace_nan_with_mean(df):
    for i in range(len(stacked_dfs_nan)):
        if i == 0:
            df[i] = df[i+1]
        if i == len(stacked_dfs_nan)-1:
            df[i] = df[i-1]
        if stacked_dfs_nan[i]:
            df[i] = (df[i-1] + df[i+1])/2
    return df

stacked_dfs = replace_nan_with_mean(stacked_dfs)
# Prepare stacke_dfs for usage in an LSTM, by creating an input window of 6 * 24 hours, and an output of 24 hours
window_size_input = 24*6
window_size_output = 24
def df_to_X_y(df, window_size_input=window_size_input, window_size_output=window_size_output):
    df_as_np = df.to_numpy()
    X = []
    y = []
    for i in range(len(df_as_np) - window_size_input-window_size_output):
        X.append([[a] for a in df_as_np[i:i+window_size_input]])
        y.append(df_as_np[i+window_size_input:i+window_size_input+window_size_output])
    print(np.array(X).shape, np.array(y).shape)
    return np.array(X), np.array(y)

X,y = df_to_X_y(stacked_dfs)

X_train, y_train = X[:15000], y[:15000]
X_val, y_val = X[15000:18000], y[15000:18000]
X_test, y_test = X[18000:], y[18000:]

# Load the model
model = load_model(r"model1/")

train_pred = model.predict(X_train)
# print shape of train_pred and y_train to see if they match

train_results = [pd.DataFrame(data={'Train Predictions': train_pred[i*24], 'True Values': y_train[i*24]}) for i in range(9)]

# test predictions
test_pred = model.predict(X_test)
test_results = [pd.DataFrame(data={'Test Predictions': test_pred[i*24], 'True Values': y_test[i*24]}) for i in range(9)]
print(train_pred.shape, y_train.shape, test_pred.shape, y_test.shape)
# create 3 x 3 subplots of the first 9 training and test predictions
fig, axs = plt.subplots(3,3, figsize=(16,9))
print(axs)
for i in range(3):
    for j in range(3):
        ax = axs[i,j]
        ax.plot(train_results[i*3 + j]['Train Predictions'], label='Train Predictions')
        ax.plot(train_results[i*3 + j]['True Values'], label='Train True Values')
        ax.plot(test_results[i*3 + j]['Test Predictions'], label='Test Predictions')
        ax.plot(test_results[i*3 + j]['True Values'], label='Test True Values')
        ax.legend()
plt.show()