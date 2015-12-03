# need to pass it a file, where data starts,  path to write

# things to import
import sys
import pandas
import scipy
import numpy
from scipy import stats
from scipy.stats import t
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
# arguments being passed
path_of_file=sys.argv[1]
last_metadata_column=int(sys.argv[2])
old_name=str(sys.argv[3])
path_to_write=sys.argv[4]

# spearman p calc based on two tailed t-test
def spearmanp(r,n):
    tstat=r*numpy.sqrt((n-2)/(1-r**2))
    return t.cdf(-abs(tstat),n-2)*2

# read in the data
df=pandas.read_table(path_of_file,index_col=False)
df.rename(columns={old_name : 'splitter'}, inplace=True)
#df=pandas.read_table("~/Desktop/temp_df.txt",index_col=False)
splitters=numpy.unique(df.splitter)
df_final=pandas.DataFrame(columns=["otus","variable","value","p.value","splitter"])

for numb in range(0,(len(splitters))):
    df_sub=df.loc[(df.splitter==splitters[numb])]
    df_data_only=df.drop(df_sub.columns[[range(0,last_metadata_column)]],axis=1)
    df_corr_matrix=df_data_only.corr(method="spearman")
    df_corr_matrix["otus"]=df_corr_matrix.index
    #melt dataframe but maintain indices now called otus
    df_melt=pandas.melt(df_corr_matrix,id_vars="otus")
    # remove NAs or NaNs which are result of non-existent otus (all 0 values)
    df_melt=df_melt[numpy.isfinite(df_melt.value)]
    df_melt['pvalue']=spearmanp(df_melt.value,df_sub.shape[0])
    df_melt['fdr'] = stats.p_adjust(FloatVector(df_melt.pvalue), method = 'fdr')
    df_melt_positive=df_melt[df_melt['value'] > 0 & df_melt['fdr'] < 0.05]

    df_melt_positive['splitter']=splitters[numb]
    df_final=df_final.append(df_melt_positive,ignore_index=True)
    df_final.to_csv((path_to_write+splitters[numb]+"final_co-occurrence_results.csv"),index=False)
#write the file
# df_final.to_csv(path_to_write,index=False)
