# co-occurrence in python
Some of my co-occurrence code with just spearman's correlations can run slowly in R.
I decided to rewrite it in python.
The script requires 4 arguments and a dataset where the first few columns are metadata and the rest are otu counts.  Below are the arguments:

1) path to data (~/data/this_is_my_data.txt)

2) last column where metadata exists (if 1 through 10 are metadata, you would give 10 as an argument)

3) name of the column that contains factors that networks will be split across (e.g. "treatments")

4) path to output and name of output (~/data/this_is_my_output.csv)

An example of what the execution looks like:
`
python ~/first_attempt_w_splitter.py ~/Desktop/temp_df.txt 7 treatments ~/Desktop/new_temp.csv
`
Note that this script requires scipy, pandas, and numpy
