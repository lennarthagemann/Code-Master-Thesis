import os
from math import sqrt
from tensorflow.keras.losses import MeanSquaredError
from statsmodels.tsa.stattools import adfuller
from statsmodels.tsa.arima.model import ARIMA
from statsmodels.graphics.tsaplots import plot_predict
from statsmodels.graphics.api import qqplot
import statsmodels.api as sm
from pmdarima import auto_arima
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from datetime import datetime, timedelta

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
print(stacked_dfs)


def plot_all_prices(prices_df):
    ax = prices_df.plot()
    ax.set_title("Historic Price Data", fontweight = "bold", fontsize = 16)
    ax.set_xlabel("Date", fontweight = "bold", fontsize = 14)
    ax.set_ylabel("Price [€/MWh]", fontweight = "bold", fontsize = 14)
    plt.show()
    return

def ad_test(dataset):
     dftest = adfuller(dataset, autolag = 'AIC')
     print("1. ADF : ",dftest[0])
     print("2. P-Value : ", dftest[1])
     print("3. Num Of Lags : ", dftest[2])
     print("4. Num Of Observations Used For ADF Regression:",      dftest[3])
     print("5. Critical Values :")
     for key, val in dftest[4].items():
         print("\t",key, ": ", val)

fig = plt.figure(figsize=(12, 8))
ax1 = fig.add_subplot(211)
fig = sm.graphics.tsa.plot_acf(stacked_dfs.values.squeeze(), lags=40, ax=ax1)
ax2 = fig.add_subplot(212)
fig = sm.graphics.tsa.plot_pacf(stacked_dfs, lags=40, ax=ax2)
plt.show()


train = stacked_dfs.iloc[:-1000]
test = stacked_dfs.iloc[-1000:]


arma_mod315 = ARIMA(train, order=(3,1,5)).fit()
print(arma_mod315.aic, arma_mod315.bic, arma_mod315.hqic)
print(arma_mod315.params)
print(sm.stats.durbin_watson(arma_mod315.resid.values))

fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111)
ax = arma_mod315.resid.plot(ax=ax)
resid = arma_mod315.resid
stats.normaltest(resid)


fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111)
fig = qqplot(resid, line="q", ax=ax, fit=True)

fig = plt.figure(figsize=(12, 8))
ax1 = fig.add_subplot(211)
fig = sm.graphics.tsa.plot_acf(resid.values.squeeze(), lags=40, ax=ax1)
ax2 = fig.add_subplot(212)
fig = sm.graphics.tsa.plot_pacf(resid, lags=40, ax=ax2)

fig, ax = plt.subplots(figsize=(10, 8))
fig = plot_predict(arma_mod315, start=datetime.strptime("2023-02-25", "%Y-%m-%d") , end=datetime.strptime("2023-04-25", "%Y-%m-%d") , ax=ax)
plt.show()

# ad_test(stacked_dfs) # Prices are appropriate for analysis with ARIMA model
# stepwise_fit = auto_arima(stacked_dfs, trace=True)  #Best Fit:Arima(3,1,5)(0,0,0)[0]

sim = arma_mod315.simulate(24, anchor='end', repetitions=100)
fig, ax = plt.subplots(figsize=(15, 4))
train.plot(color='k', ax=ax)
sim.plot(ax=ax, color='C0', alpha=0.2, legend=False)
plt.show()