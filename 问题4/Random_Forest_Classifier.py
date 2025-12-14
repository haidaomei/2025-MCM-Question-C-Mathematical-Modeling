import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score

plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False

df = pd.read_csv('女胎检测数据.csv', header=0, encoding="gbk")

feature_columns = [10, 16, 17, 18, 19, 22, 23, 24, 25]
X = df.iloc[0:605, feature_columns]

y = df.iloc[0:605, 27]

y = y.notna().astype(int)

print(f"正常样本数(y=0):{(y==0).sum()}\n")
print(f"异常样本数(y=1):{(y==1).sum()}\n")

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

print(f"训练集大小:{X_train.shape[0]}\n")
print(f"测试集大小:{X_test.shape[0]}\n")

model = RandomForestClassifier(n_estimators=1000)

model.fit(X_train, y_train)

y_pred = model.predict(X_test)

accuracy = accuracy_score(y_test, y_pred)
print(f"测试集准确率:{accuracy}\n")

feature_names = [f"列{col}" for col in feature_columns]
print(f"特征重要性\n")
for name, importance in zip(feature_names, model.feature_importances_):
    print(f"{name}:{importance}\n")

#for i in range(len(X_test)):
#    true_label = "正常" if y_test.iloc[i] == 0 else "异常"
#    pred_label = "正常" if y_pred[i] == 0 else "异常"
#    print(f"样本{i+1}:真实={true_label},预测={pred_label}\n")

correct_predictions = (y_pred == y_test).astype(int)

colors = ['red' if cp == 0 else 'green' for cp in correct_predictions]

plt.scatter(range(len(y_test)), correct_predictions, c=colors)

plt.xlabel('测试样本序号')
plt.ylabel('预测结果')
plt.title(f'随机森林模型预测结果(准确率:{accuracy})')
plt.yticks([0, 1], ['预测错误', '预测正确'])

from matplotlib.patches import Patch

legend_elements = \
[
    Patch(facecolor='green', label='预测正确'),
    Patch(facecolor='red', label='预测错误')
]
plt.legend(handles=legend_elements, loc='upper right')

plt.savefig('随机森林模型预测结果.png')

plt.show()