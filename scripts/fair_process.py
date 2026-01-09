import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression

data = {
    "N_monomer": [1,2,3,4,8],
    "Delta_G":   [-2.85060894,-4.615932618,-6.376827417,-8.138869297,-15.17875589],
}

df = pd.DataFrame(data)

X = df[["N_monomer"]].values
y = df["Delta_G"].values

model = LinearRegression()
model.fit(X, y)

a = model.coef_[0]
b = model.intercept_

print("拟合得到：ΔG_form(N) ≈ a * N + b")
print("a (per monomer) =", a)
print("b (end effect)  =", b)

df["DeltaG_pred"] = model.predict(X)
df["residual"] = df["Delta_G"] - df["DeltaG_pred"]

print(df)
