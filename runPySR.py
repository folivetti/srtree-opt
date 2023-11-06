import pysr
import pandas as pd 

df = pd.read_csv("SeoulBikeDataNZ.csv")

x = df[["Temperature", "Humidity", "WindSpeed", "Visibility", "SolarRadiation", "date_dec_yearly", "FunctioningDay"]].values
y = df.RentedBikeCount.values

myloss = """
function eval_loss(tree, dataset::Dataset{T,L}, options)::L where {T,L}
    prediction, flag = eval_tree_array(tree, dataset.X, options)
    if !flag
        return L(Inf)
    end
    return sum(exp.(exp.(prediction) .- dataset.y .* prediction .- dataset.y))
end
"""

model = pysr.PySRRegressor(loss="myloss(yhat, ys) = sum(exp(yhat) - ys * yhat + ys*log(ys))", unary_operators=["sqrt", "log1p"], nested_constraints={"sqrt" : {"sqrt" : 0, "log1p" : 0}, "log1p" : {"sqrt" : 0, "log1p" : 0}}, population_size=100, batching=False)
#model = pysr.PySRRegressor(full_objective=myloss, population_size=500)
model.fit(x,y)

#model = pysr.PySRRegressor(loss="myloss(yhat, ys) = exp(exp(yhat) - ys * yhat - ys)", unary_operators=["log1p"])
#model.fit(x,y)
