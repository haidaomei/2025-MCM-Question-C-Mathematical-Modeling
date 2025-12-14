import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False

df = pd.read_csv('男胎检测数据.csv', encoding='gbk')

data = df.iloc[1:1082]

bmi = data.iloc[:, 11].astype(float)
y_concentration = data.iloc[:, 22].astype(float)
weeks = data.iloc[:, 10].astype(float)

bmi_groups = \
[
    (20, 25),
    (25, 42),
    (42, 50)
]

print("朴素分组分析\n")
print("各BMI分组中Y染色体浓度在0.04-0.1范围内的检测孕周平均值:\n")

for i, (bmi_low, bmi_high) in enumerate(bmi_groups, 1):
    bmi_mask = (bmi >= bmi_low) & (bmi <= bmi_high)

    y_mask = (y_concentration >= 0.04) & (y_concentration <= 0.1)

    mask = bmi_mask & y_mask

    filtered_weeks = weeks[mask]

    if len(filtered_weeks) > 0:
        mean_week = filtered_weeks.mean()
        print(f"组{i}(BMI {bmi_low}-{bmi_high}):{mean_week}周\n")

X = bmi.values.reshape(-1, 1)

kmeans = KMeans(n_clusters=6)
clusters = kmeans.fit_predict(X)

data_copy = data.copy()
data_copy['BMI'] = bmi
data_copy['Y染色体浓度'] = y_concentration
data_copy['检测孕周'] = weeks
data_copy['cluster'] = clusters

import matplotlib.pyplot as plt
for cluster_id in range(6):
    cluster_data = data_copy[data_copy['cluster'] == cluster_id]
    plt.scatter(cluster_data['BMI'], cluster_data['Y染色体浓度'], label=f'聚类{cluster_id+1}')

plt.xlabel('BMI')
plt.ylabel('Y染色体浓度')
plt.title('BMI聚类分布')
plt.legend()
plt.savefig('单元K-Means.png')
plt.show()

print("KMeans聚类分析\n")
print("特征:孕妇BMI\n")
print("各聚类中Y染色体浓度在0.04-0.1范围内的检测孕周平均值:\n")
for cluster_id in range(6):
    cluster_data = data_copy[data_copy['cluster'] == cluster_id]

    filtered_data = cluster_data[(cluster_data['Y染色体浓度'] >= 0.04) & (cluster_data['Y染色体浓度'] <= 0.1)]

    if len(filtered_data) > 0:
        mean_week = filtered_data['检测孕周'].mean()
        print(f"聚类{cluster_id+1}:{mean_week}周\n")

print("各聚类BMI区间:\n")
for cluster_id in range(6):
    cluster_data = data_copy[data_copy['cluster'] == cluster_id]
    min_bmi = cluster_data['BMI'].min()
    max_bmi = cluster_data['BMI'].max()
    print(f"聚类{cluster_id+1}:BMI范围[{min_bmi},{max_bmi}]\n")