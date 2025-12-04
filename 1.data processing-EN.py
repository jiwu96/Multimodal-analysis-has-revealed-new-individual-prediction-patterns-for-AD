import os
import re
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multitest import multipletests
from combat.pycombat import pycombat

os.chdir("D:\\AD自测数据\\GEO-中国的")

df = pd.read_csv("GSE239454.txt", sep="\t", header=0)
df.iloc[:, 0] = df.iloc[:, 0].str.replace(r"ENSG\d+\|", "", regex=True)
df.to_csv("output.txt", sep="\t", index=False, quote=False)

rt = pd.read_csv("output.txt", sep="\t", index_col=0)
exp = rt.values
data = pd.DataFrame(exp, index=rt.index, columns=rt.columns)

qx = np.quantile(data.values.flatten(), [0, 0.25, 0.5, 0.75, 0.99, 1.0])
LogC = (qx[4] > 100) or ((qx[5] - qx[0] > 50) and qx[1] > 0)
if LogC:
    data[data < 0] = 0
    data = np.log2(data + 1)

scaler = StandardScaler()
data_normalized = pd.DataFrame(
    scaler.fit_transform(data.T).T,
    index=data.index,
    columns=data.columns
)

con = pd.read_csv("s1.txt", sep="\t", header=None)
treat = pd.read_csv("s2.txt", sep="\t", header=None)
conData = data_normalized[con[0]]
treatData = data_normalized[treat[0]]
data_combined = pd.concat([conData, treatData], axis=1)
conNum = conData.shape[1]
treatNum = treatData.shape[1]

Type = ["Control"] * conNum + ["Treat"] * treatNum
outData = pd.concat([pd.Series([f"{col}_{typ}" for col, typ in zip(data_combined.columns, Type)]), 
                     data_combined.reset_index(drop=True)])
outData.to_csv("GSE239454.normalize.txt", sep="\t", index=False, header=False)

sampleInfo = pd.DataFrame({
    'sampleName': data_combined.columns,
    'group': Type
})

trainSamples, testSamples = train_test_split(
    sampleInfo['sampleName'], 
    test_size=0.3, 
    stratify=sampleInfo['group'],
    random_state=123
)

trainExpr = data_combined[trainSamples]
testExpr = data_combined[testSamples]

trainOut = pd.concat([pd.Series([f"{col}_{sampleInfo[sampleInfo['sampleName'] == col]['group'].values[0]}" 
                                for col in trainExpr.columns]), 
                      trainExpr.reset_index(drop=True)])
testOut = pd.concat([pd.Series([f"{col}_{sampleInfo[sampleInfo['sampleName'] == col]['group'].values[0]}" 
                               for col in testExpr.columns]), 
                     testExpr.reset_index(drop=True)])

trainOut.to_csv("train.70.txt", sep="\t", index=False, header=False)
testOut.to_csv("test.30.txt", sep="\t", index=False, header=False)

print(f"拆分完成：\n训练集样本数：{len(trainSamples)}\n验证集样本数：{len(testSamples)}")

os.chdir("D:\\AD自测数据\\所有数据合并\\1")
files = [f for f in os.listdir() if f.endswith(".txt") and f not in ["merge.preNorm.txt", "merge.normalize.txt"]]
geneList = {}

for file in files:
    rt = pd.read_csv(file, sep="\t", index_col=0)
    geneNames = rt.index.tolist()
    uniqGene = list(set(geneNames))
    header = re.split(r"\.|-", file)
    geneList[header[0]] = uniqGene

interGenes = set.intersection(*[set(genes) for genes in geneList.values()])

allTab = pd.DataFrame()
batchType = []
for i, file in enumerate(files):
    rt = pd.read_csv(file, sep="\t", index_col=0)
    exp = rt.values
    data = pd.DataFrame(exp, index=rt.index, columns=rt.columns)
    data = data.groupby(data.index).mean()
    header = re.split(r"\.|-", file)
    data.columns = [f"{header[0]}_{col}" for col in data.columns]
    
    if i == 0:
        allTab = data.loc[interGenes]
    else:
        allTab = pd.concat([allTab, data.loc[interGenes]], axis=1)
    
    batchType.extend([i] * data.shape[1])

outTab = pd.concat([pd.Series(allTab.columns), allTab.reset_index(drop=True)])
outTab.to_csv("merge.preNorm.txt", sep="\t", index=False, header=False)

allTab_corrected = pycombat(allTab, batchType)
outTab_corrected = pd.concat([pd.Series(allTab_corrected.columns), allTab_corrected.reset_index(drop=True)])
outTab_corrected.to_csv("merge.normalize.txt", sep="\t", index=False, header=False)

def bioPCA(inputFile, outFile, titleName):
    rt = pd.read_csv(inputFile, sep="\t", index_col=0)
    data = rt.T
    Project = data.index.str.extract(r"(.*?)\_.*")[0]
    
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(data)
    PCA_df = pd.DataFrame({'PC1': pca_result[:, 0], 'PC2': pca_result[:, 1], 'Type': Project})
    
    plt.figure(figsize=(5.5, 4.25))
    sns.scatterplot(x='PC1', y='PC2', hue='Type', style='Type', data=PCA_df)
    plt.title(titleName)
    plt.savefig(outFile)
    plt.close()

bioPCA("merge.preNorm.txt", "PCA.preNorm.pdf", "Before batch correction")
bioPCA("merge.normalize.txt", "PCA.normalzie.pdf", "After batch correction")

os.chdir("D:\\AD自测数据\\所有数据合并\\1\\T_cell")
expFile = "merge.normalize.txt"
geneFile = "gene.txt"

rt = pd.read_csv(expFile, sep="\t", index_col=0)
data = rt.T
Type = data.index.str.extract(r".*\\_(.*?)")[0]
treatData = data[Type == "Treat"]
rt = pd.concat([data, Type], axis=1)

data_melted = rt.melt(id_vars=["Type"], var_name="Gene", value_name="Expression")

plt.figure(figsize=(45, 4.5))
sns.boxplot(x="Gene", y="Expression", hue="Type", data=data_melted)
plt.xticks(rotation=45)
plt.savefig("boxplot.pdf")
plt.close()

p_values = []
genes = data_melted['Gene'].unique()
for gene in genes:
    control = data_melted[(data_melted['Gene'] == gene) & (data_melted['Type'] == 'Control')]['Expression']
    treat = data_melted[(data_melted['Gene'] == gene) & (data_melted['Type'] == 'Treat')]['Expression']
    p_value = stats.mannwhitneyu(control, treat).pvalue
    p_values.append(p_value)

reject, pvals_corrected, _, _ = multipletests(p_values, method='fdr_bh')
significant_genes = genes[reject]

pd.Series(significant_genes).to_csv("significant_genes.txt", index=False, header=False)

geneFile = "significant_genes.txt"
files = [f for f in os.listdir() if f.endswith("normalize.txt")]
geneList = {}

for file in files:
    rt = pd.read_csv(file, sep="\t", index_col=0)
    geneNames = rt.index.tolist()
    uniqGene = list(set(geneNames))
    header = re.split(r"\.|-", file)
    geneList[header[0]] = uniqGene

interGenes = set.intersection(*[set(genes) for genes in geneList.values()])

allTab = pd.DataFrame()
batchType = []
for i, file in enumerate(files):
    rt = pd.read_csv(file, sep="\t", index_col=0)
    data = rt.groupby(rt.index).mean()
    header = re.split(r"\.|-", file)
    data.columns = [f"{header[0]}_{col}" for col in data.columns]
    
    if i == 0:
        allTab = data.loc[interGenes]
    else:
        allTab = pd.concat([allTab, data.loc[interGenes]], axis=1)
    
    batchType.extend([i] * data.shape[1])

svaTab = pycombat(allTab, batchType)
geneRT = pd.read_csv(geneFile, header=None)
geneTab = svaTab.loc[intersection(svaTab.index, geneRT[0])].T

train = geneTab.index.str.startswith("merge")
trainExp = geneTab[train]
testExp = geneTab[~train]
trainExp.index = trainExp.index.str.replace("merge_", "Train.")
testExp.index = testExp.index.str.replace(".*_", "", regex=True)

trainType = trainExp.index.str.extract(r".*\\_(.*)")[0]
testType = testExp.index.str.extract(r".*\\_(.*)")[0]
trainType = np.where(trainType == "Control", 0, 1)
testType = np.where(testType == "Control", 0, 1)

trainExp['Type'] = trainType
testExp['Type'] = testType

trainExp.to_csv("data.train.txt", sep="\t", index=False)
testExp.to_csv("data.test.txt", sep="\t", index=False)