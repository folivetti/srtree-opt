from pyoperon.sklearn import SymbolicRegressor
import pandas as pd 

#df = pd.read_csv("SeoulBikeDataSummer.csv")
#df = pd.read_csv("wine.csv")
df = pd.read_csv("apple1.csv").dropna()
x = df[['regprc', 'ecoprc', 'inseason', 'hhsize', 'faminc', 'age']].values
y = df.reglbs.values

#x = df[["Temperature", "Humidity", "WindSpeed", "Visibility", "SolarRadiation", "date_dec_yearly", "FunctioningDay"]].values
#x = df[["alcohol", "deaths", "liver"]].values
#y = df.RentedBikeCount.values
#y = df.heart.values


model = SymbolicRegressor(optimizer="lbfgs", optimizer_likelihood="poisson", optimizer_likelihood_loginput=False, max_length=50, objectives=['r2','length'], allowed_symbols="add,mul,constant,variable,sqrt,log")
model.fit(x,y)
print(model.get_model_string(model.model_))
print(model)

