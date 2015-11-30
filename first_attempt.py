import pandas
import scipy
import numpy
from scipy import stats
from scipy.stats import t
def spearmanp(r,n):
    tstat=r*numpy.sqrt((n-2)/(1-r**2))
    return t.cdf(-abs(tstat),n-2)*2






df=pandas.read_table("~/Desktop/pitfoam_analysis/merged_pitfoam_table.txt",index_col=False)
df_sub=df.loc[(df.splitter=="foaming_untreated_consistently foaming")]
df_data_only=df_sub.drop(df_sub.columns[[range(0,7)]],axis=1)
df_corr_matrix=df_data_only.corr(method="spearman")
