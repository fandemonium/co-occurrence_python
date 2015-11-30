# need to pass it a file, where data starts,  path to write

# things to import
import sys
import pandas
import scipy
import numpy
from scipy import stats
from scipy.stats import t

# arguments being passed
path_of_file=sys.argv[1]
last_metadata_column=int(sys.argv[2])
path_to_write=sys.argv[3]

# spearman p calc based on two tailed t-test
def spearmanp(r,n):
    tstat=r*numpy.sqrt((n-2)/(1-r**2))
    return t.cdf(-abs(tstat),n-2)*2

# read in the data
df=pandas.read_table(path_of_file,index_col=False)
# remove metadata columns
df_data_only=df.drop(df.columns[[range(0,last_metadata_column)]],axis=1)
#make correlation matrix
df_corr_matrix=df_data_only.corr(method="spearman")
#make column based on rows (called indexes in python)
df_corr_matrix["otus"]=df_corr_matrix.index
#melt dataframe but maintain indices now called otus
df_melt=pandas.melt(df_corr_matrix,id_vars="otus")
# remove NAs or NaNs which are result of non-existent otus (all 0 values)
df_melt=df_melt[numpy.isfinite(df_melt.value)]
df_melt['p.value']=spearmanp(df_melt.value,df_sub.shape[0])
#write the file
df_melt.to_csv(path_to_write,index=False)
