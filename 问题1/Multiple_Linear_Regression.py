import pandas as pd
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False

df = pd.read_csv('男胎检测数据.csv', encoding='gbk')

variables_data = {}
col_indices = \
{
    'X染色体浓度': 23,
    'Y染色体浓度': 22,
    '18号染色体的Z值': 18,
    '在参考基因组上比对的比例': 13,
    '孕妇BMI': 11,
    '年龄': 2,
    'Y染色体的Z值': 21,
    '重复读段的比例': 14
}

for var_name, col_idx in col_indices.items():
    data = []
    for i in range(1, 1082):
        data.append(float(df.iloc[i, col_idx]))
    variables_data[var_name] = np.array(data)

Y = variables_data['Y染色体浓度']

X_vars = ['X染色体浓度', '18号染色体的Z值', '在参考基因组上比对的比例', '孕妇BMI', '年龄', 'Y染色体的Z值', '重复读段的比例']

print("多元线性回归分析:X染色体浓度 18号染色体的Z值 在参考基因组上比对的比例 孕妇BMI 年龄 Y染色体的Z值 重复读段的比例 Y染色体浓度")

X = np.column_stack([variables_data[var] for var in X_vars])

X_with_const = sm.add_constant(X)

model = sm.OLS(Y, X_with_const).fit()

print(model.summary())

Y_pred = model.predict(X_with_const)
residuals = Y - Y_pred

SS_res = np.sum(residuals**2)
R2 = model.rsquared

print("回归方程:")
equation = f"Y染色体浓度={model.params[0]:}"
for i, var_name in enumerate(X_vars, 1):
    coef_sign = "+" if model.params[i] >= 0 else "-"
    equation += f"{coef_sign}{abs(model.params[i])}×{var_name}"
print(equation)

print(f"R²={R2}")
print(f"Q={SS_res}")

plt.scatter(Y_pred, Y)
plt.plot([Y.min(), Y.max()], [Y.min(), Y.max()])
plt.xlabel('预测值')
plt.ylabel('实际值')
plt.title('实际值 预测值')
plt.savefig('多元线性回归分析 实际值 预测值.png')
plt.show()

plt.scatter(Y_pred, residuals)
plt.axhline(y=0)
plt.xlabel('预测值')
plt.ylabel('残差')
plt.title('残差 预测值')
plt.savefig('多元线性回归分析 残差 预测值.png')
plt.show()

plt.hist(residuals, bins=100)
plt.xlabel('残差')
plt.ylabel('频数')
plt.title('残差分布')
plt.savefig('多元线性回归分析 残差分布.png')
plt.show()

from scipy import stats

shapiro_stat, shapiro_p = stats.shapiro(residuals)
print(f"Shapiro-Wilk检验统计量:p值:{shapiro_p}")