from pyoperon.sklearn import SymbolicRegressor
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


model = SymbolicRegressor(optimizer="lbfgs", optimizer_likelihood="poisson", optimizer_likelihood_loginput=True, max_length=20, objectives=['r2','length'])
model.fit(x,y)
print(model.get_model_string(model.model_))
print(model)

