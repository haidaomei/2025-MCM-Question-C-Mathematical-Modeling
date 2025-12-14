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

print("单元线性回归分析:X染色体浓度 Y染色体浓度")

X1 = sm.add_constant(X_chr)
model1 = sm.OLS(Y_chr, X1).fit()
print(model1.summary())

Y_pred = model1.predict(X1)
residuals = Y_chr - Y_pred

Q = np.sum(residuals**2)

print(f"回归方程:y={model1.params[0]:}+{model1.params[1]:}x")
print(f"R²={model1.rsquared:}")
print(f"Q={Q}")

plt.scatter(X_chr, Y_chr)
x_line = np.linspace(X_chr.min(), X_chr.max())
y_line1 = model1.params[0] + model1.params[1] * x_line
plt.plot(x_line, y_line1)
plt.xlabel('X染色体浓度')
plt.ylabel('Y染色体浓度')
plt.title('X染色体浓度 Y染色体浓度')
plt.savefig('单元线性回归分析 X染色体浓度 Y染色体浓度.png')
plt.show()

plt.scatter(Y_pred, residuals)
plt.axhline(y=0)
plt.xlabel('预测值')
plt.ylabel('残差')
plt.title('残差 预测值')
plt.savefig('单元线性回归分析 残差 预测值.png')
plt.show()

plt.hist(residuals, bins=100)
plt.xlabel('残差')
plt.ylabel('频数')
plt.title('残差分布')
plt.savefig('单元线性回归分析 残差分布.png')
plt.show()

from scipy import stats

shapiro_stat, shapiro_p = stats.shapiro(residuals)
print(f"Shapiro-Wilk检验p值:{shapiro_p}")