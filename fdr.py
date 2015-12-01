import sys
import pandas
import scipy
import numpy
from scipy import stats
from scipy.stats import t
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector

path_of_file=sys.argv[1]
old_name=str(sys.argv[2])
cutoff=float(sys.argv[3])
path_to_write=sys.argv[4]


stats=importr('stats')
df=pandas.read_csv(path_of_file,index_col=False)
df.rename(columns={old_name: 'pvalue'}, inplace=True)
df['fdr'] = stats.p_adjust(FloatVector(df.pvalue), method = 'fdr')
df_sub=df[df['fdr'] < cutoff]
df_sub.to_csv(path_to_write,index=False)
