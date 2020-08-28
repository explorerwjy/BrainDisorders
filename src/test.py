from utils import *
from IPython.display import display, HTML
import os


a, b = 3.23, 1.12
LGD_Dat = LoadLGDDat("./NMD_model/BRNAMY.obs")
post = ComputeLikelihoodNMDPara(LGD_Dat, a, b, math.log(scipy.stats.halfcauchy.pdf(1/math.sqrt(a+b))))
print(post)
