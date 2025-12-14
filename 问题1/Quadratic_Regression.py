import pandas as pd
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False

df = pd.read_csv('男胎检测数据.csv', encoding='gbk')

X_chr_data = []
Y_chr_data = []

for i in range(1, 1082):
    X_chr_data.append(float(df.iloc[i, 23]))
    Y_chr_data.append(float(df.iloc[i, 22]))

X_chr = np.array(X_chr_data)
Y_chr = np.array(Y_chr_data)

print("单元二次回归分析:X染色体浓度 Y染色体浓度")

X_poly = np.column_stack([np.ones(len(X_chr)), X_chr, X_chr**2])

model_poly = sm.OLS(Y_chr, X_poly).fit()

print(model_poly.summary())

Y_pred = model_poly.predict(X_poly)
residuals = Y_chr - Y_pred

SS_res = np.sum(residuals**2)
R2 = model_poly.rsquared

print("回归方程:")
equation = f"y={model_poly.params[0]:}"

var_names_poly = ['x', 'x²']

for i, var_name in enumerate(var_names_poly, 1):
    coef_sign = "+" if model_poly.params[i] >= 0 else "-"
    equation += f"{coef_sign}{abs(model_poly.params[i])}{var_name}"

print(equation)

print(f"R²={R2}")
print(f"Q={SS_res}")

plt.scatter(X_chr, Y_chr)
x_line = np.linspace(X_chr.min(), X_chr.max(), 100)
x_line_design = np.column_stack([np.ones(100), x_line, x_line**2])
y_line_poly = model_poly.predict(x_line_design)
plt.plot(x_line, y_line_poly)
plt.xlabel('X染色体浓度')
plt.ylabel('Y染色体浓度')
plt.title('X染色体浓度 Y染色体浓度')
plt.savefig('单元二次回归分析 X染色体浓度 Y染色体浓度.png')
plt.show()

plt.scatter(Y_pred, residuals)
plt.axhline(y=0)
plt.xlabel('预测值')
plt.ylabel('残差')
plt.title('残差 预测值')
plt.savefig('单元二次回归分析 残差 预测值.png')
plt.show()

plt.hist(residuals, bins=100)
plt.xlabel('残差')
plt.ylabel('频数')
plt.title('残差分布')
plt.savefig('单元二次回归分析 残差分布.png')
plt.show()

from scipy import stats

shapiro_stat, shapiro_p = stats.shapiro(residuals)
print(f"Shapiro-Wilk检验p值:{shapiro_p}")