import pandas as pd
from sklearn.feature_selection import f_classif
from statsmodels.stats.multitest import fdrcorrection


data = pd.read_csv("trans_for_ftest.csv")

y = data["status"]
x = data.drop("status", axis=1)

sel = f_classif(x, y)

p_values = pd.Series(sel[1])

p_val = fdrcorrection(p_values, 0.05) 
p_val_adjusted = p_val[1]

p_val_adjusted = pd.Series(p_val_adjusted)

p_val_adjusted.index = x.columns 

print(p_val_adjusted[p_val_adjusted<0.01])


