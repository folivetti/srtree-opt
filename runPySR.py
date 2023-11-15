import pysr
import pandas as pd 

#df = pd.read_csv("SeoulBikeDataSummer.csv")
#df = pd.read_csv("wine.csv")
df = pd.read_csv("fertil2.csv").dropna()
df = df[df.children != 0]

#x = df[["Temperature", "Humidity", "WindSpeed", "Visibility", "SolarRadiation", "date_dec_yearly", "FunctioningDay"]].values
#x = df[["alcohol", "deaths", "liver"]].values
x = df[['educ','age','evermarr','urban','electric','tv']].values
#y = df.RentedBikeCount.values
#y = df.heart.values
y = df.children.values


model = pysr.PySRRegressor(loss="myloss(yhat, ys) = sum(exp(yhat) - ys * yhat + ys*log(ys))", unary_operators=["sqrt", "log1p"], nested_constraints={"sqrt" : {"sqrt" : 0, "log1p" : 0}, "log1p" : {"sqrt" : 0, "log1p" : 0}}, population_size=100, maxsize=25, batching=False)
#model = pysr.PySRRegressor(loss="myloss(yhat, ys) = sum(exp(yhat) - ys * yhat + ys*log(ys))", population_size=500, maxsize=50, batching=False, nested_constraints={"/" : {"/" : 1}})
model.fit(x,y)

#model = pysr.PySRRegressor(loss="myloss(yhat, ys) = exp(exp(yhat) - ys * yhat - ys)", unary_operators=["log1p"])
#model.fit(x,y)
