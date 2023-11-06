# coding: utf-8
from pyoperon.sklearn import SymbolicRegressor
import pandas as pd
import numpy as np

df = pd.read_csv("SeoulBikeData.csv")
X = df[['Date_dec', 'date_dec_yearly', 'date_dec_weekly', 
       'Hour', 'Temperature', 'Humidity', 'WindSpeed', 'Visibility',
       'Visibility_2000_or_more', 'DewPointTemperature', 'SolarRadiation',
       'Rainfall', 'Snowfall', 'Seasons_Winter', 'Seasons_Summer',
       'Seasons_Autumn', 'Seasons_Spring', 'Holiday', 'FunctioningDay']].values
y = df.RentedBikeCount.values
reg = SymbolicRegressor(max_length=25, optimizer="lbfgs", optimizer_likelihood="poisson", optimizer_likelihood_loginput=False, objectives=["mse", "length"], allowed_symbols="add,sub,mul,div,log,constant,variable", random_state=42)
reg.fit(X, y)
#print(reg.get_model_string(reg.model_, precision=4))
print(reg.pareto_front_[0])
yhat = np.exp(329.952362 + (16.408461 * (1.772323 * X[:,4])))
print(yhat)
