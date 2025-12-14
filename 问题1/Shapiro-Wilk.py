import pandas as pd
from scipy import stats

df = pd.read_csv('男胎检测数据.csv',encoding='gbk')

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

print("Shapiro-Wilk正态性检验结果:\n")

for col_idx, col_name in columns_to_test.items():
    data = df.iloc[1:1082, col_idx]

    stat, p_value = stats.shapiro(data)

    alpha = 0.05
    if p_value > alpha:
        result = "符合正态分布"
    else:
        result = "不符合正态分布"

    print(f"{col_name},p值:{p_value},{result}\n")