from utils import *

ExonExp = pd.read_csv("../data/expression/brainspan/exons_matrix/expression_matrix.csv", header=None, index_col=0)
df_n = quantileNormalize(ExonExp)
df_n = df_n.round(6)
df_n.to_csv("/Users/jiayao/Work/BrainDisorders/data/expression/brainspan/exons_matrix/qn_exons_matrix.csv", sep=",", header=None)
