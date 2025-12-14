import pandas as pd
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['Microsoft YaHei']
plt.rcParams['axes.unicode_minus'] = False

df = pd.read_csv('男胎检测数据.csv', encoding='gbk')

columns_to_test = \
{
    2: "年龄",
    3: "身高",
    4: "体重",
    10: "检测孕周",
    11: "孕妇BMI",
    12: "原始读段数",
    13: "在参考基因组上比对的比例",
    14: "重复读段的比例",
    15: "唯一比对的读段数",
    16: "GC含量",
    17: "13号染色体的Z值",
    18: "18号染色体的Z值",
    19: "21号染色体的Z值",
    20: "X染色体的Z值",
    21: "Y染色体的Z值",
    22: "Y染色体浓度",
    23: "X染色体浓度",
    24: "13号染色体的GC含量",
    25: "18号染色体的GC含量",
    26: "21号染色体的GC含量",
    27: "被过滤掉读段数的比例"
}

y_concentration = []
for i in range(1, 1082):
    val = df.iloc[i, 22]
    y_concentration.append(float(val))

y_array = np.array(y_concentration)

print("各变量与Y染色体浓度的斯佩尔曼相关性分析:\n")

results = []

for col_idx, col_name in columns_to_test.items():
    if col_idx == 22:
        continue

    data = []
    for i in range(1, 1082):
        val = df.iloc[i, col_idx]
        data.append(float(val))

    x_array = np.array(data)

    valid_pairs = []
    for x, y in zip(x_array, y_array):
        if not np.isnan(x) and not np.isnan(y):
            valid_pairs.append((x, y))

    x_vals = [p[0] for p in valid_pairs]
    y_vals = [p[1] for p in valid_pairs]

    rho, p_value = stats.spearmanr(x_vals, y_vals)

    if p_value < 0.001:
        sig = "3"
    elif p_value < 0.01:
        sig = "2"
    elif p_value < 0.05:
        sig = "1"
    else:
        sig = "0"

    results.append((col_name, rho, p_value, len(valid_pairs), sig))

results.sort(key=lambda x: x[2])

for col_name, rho, p_value, n, sig in results:
    print(f"{col_name},p值:{p_value},{sig}\n")

variable_names = [r[0] for r in results]
p_values = [r[2] for r in results]

y_pos = range(len(variable_names))
plt.barh(y_pos, p_values)

plt.yticks(y_pos, variable_names)

plt.xlabel('p值')
plt.xscale('log')

plt.title('斯皮尔曼相关性分析p值分布')

plt.savefig('斯皮尔曼相关性分析p值分布.png')

plt.show()