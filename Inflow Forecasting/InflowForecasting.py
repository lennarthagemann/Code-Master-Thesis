import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import datetime as dt
from scipy import stats
from statsmodels.tsa.arima.model import ARIMA
from pmdarima import auto_arima
from pathlib import Path