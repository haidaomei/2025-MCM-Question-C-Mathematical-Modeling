import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

df = pd.read_csv('男胎检测数据.csv', encoding='gbk')

data = df.iloc[1:1082]

bmi = data.iloc[:, 11].astype(float)
y_concentration = data.iloc[:, 22].astype(float)
weeks = data.iloc[:, 10].astype(float)

x_concentration = data.iloc[:, 23].astype(float)
z_18 = data.iloc[:, 18].astype(float)
alignment_ratio = data.iloc[:, 13].astype(float)

features_df = pd.DataFrame({
    'X染色体浓度': x_concentration,
    '18号染色体Z值': z_18,
    '比对比例': alignment_ratio,
    'BMI': bmi,
    'Y染色体浓度': y_concentration,
    '检测孕周': weeks
})

X_multi = features_df[['X染色体浓度', '18号染色体Z值', '比对比例', 'BMI']].values

scaler = StandardScaler()
X_multi_scaled = scaler.fit_transform(X_multi)

kmeans_multi = KMeans(n_clusters=6, random_state=42, n_init=10)
clusters_multi = kmeans_multi.fit_predict(X_multi_scaled)

features_df['cluster_multi'] = clusters_multi

print("多维KMeans聚类分析\n")
print("特征:X染色体浓度、18号染色体的Z值、在参考基因组上比对的比例、孕妇BMI\n")
print("各聚类中Y染色体浓度在0.04-0.1范围内的检测孕周平均值:\n")

for cluster_id in range(6):
    cluster_data = features_df[features_df['cluster_multi'] == cluster_id]

    filtered_data = cluster_data[
        (cluster_data['Y染色体浓度'] >= 0.04) &
        (cluster_data['Y染色体浓度'] <= 0.1)
        ]

    if len(filtered_data) > 0:
        mean_week = filtered_data['检测孕周'].mean()
        print(f"聚类{cluster_id + 1}:{mean_week}周\n")

print("各聚类BMI区间:\n")
for cluster_id in range(6):
    cluster_data = features_df[features_df['cluster_multi'] == cluster_id]
    min_bmi = cluster_data['BMI'].min()
    max_bmi = cluster_data['BMI'].max()
    print(f"聚类{cluster_id + 1}:BMI范围[{min_bmi},{max_bmi}]\n")