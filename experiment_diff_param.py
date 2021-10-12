import BLAST_prob
from BLAST_prob import seq_query, BLAST

import numpy as np
import pandas as pd

SEQUENCE,SEQ_PROB, NEUC_PROB, Q_PROB = seq_query.read_sequence(demo_seq_len=1000)

# experiment for word size
test_wordsize = []
for l in [7,8,9]:
    for i in range(10):
        result = {}
        result['word_size'] = l
        start_pos, query = seq_query.generate_query(SEQUENCE,SEQ_PROB,query_len=50)
        result['start'] = start_pos
        result['end'] = start_pos+50-1
    
        result_hit, t = BLAST.BLAST(NEUC_PROB,query,word_l = l, window = 30) 
        
        result['result_start'] = result_hit.db_start
        result['result_end'] = result_hit.db_end
        result['time'] = t['total']
        test_wordsize.append(result)
df_wordsize = pd.DataFrame(test_wordsize)

# experiment for window size
test_windowsize = []
for s in [20,30,40,50]:
    for i in range(10):
        result = {}
        result['window_size'] = s
        start_pos, query = seq_query.generate_query(SEQUENCE,SEQ_PROB,query_len=70)
        result['start'] = start_pos
        result['end'] = start_pos+70-1
        
        result_hit, t = BLAST.BLAST(NEUC_PROB,query,word_l = 8, window = s)
        
        result['result_start'] = result_hit.db_start
        result['result_end'] = result_hit.db_end
        result['time'] = t['total']
        test_windowsize.append(result)
df_windowsize = pd.DataFrame(test_windowsize)

# experiment for query length
test_qlen = []
for q in [40,50,60,70]:
    for i in range(10):
        result = {}
        result['query_length'] = q
        start_pos, query = seq_query.generate_query(SEQUENCE,SEQ_PROB,query_len=q)
        result['start'] = start_pos
        result['end'] = start_pos+q-1
        
        result_hit, t = BLAST.BLAST(NEUC_PROB,query,word_l = 8, window = 30)
        
        result['result_start'] = result_hit.db_start
        result['result_end'] = result_hit.db_end
        result['time'] = t['total']
        test_qlen.append(result)
df_qlen = pd.DataFrame(test_qlen)

# calculate accuracy
for df in [df_wordsize,df_windowsize,df_qlen]:
    if (df['result_start'].all()):
        df['match_start'] = np.where(df['start'] > df['result_start'], df['start'], df['result_start'])
        df['match_end'] = np.where(df['end'] < df['result_end'], df['end'], df['result_end'])
        df['accuracy'] =  (df['match_end'] - df['match_start'] + 1) / (df['end'] - df['start'] + 1)
    else:
        df['accuracy'] = 0
    
# get the mean performance for different values of each parameters
df_wordsize = df_wordsize[['word_size', 'time', 'accuracy']].groupby(['word_size']).mean()
print(df_wordsize)
df_windowsize = df_windowsize[['window_size', 'time', 'accuracy']].groupby(['window_size']).mean()
print(df_windowsize)
df_qlen = df_qlen[['query_length', 'time', 'accuracy']].groupby(['query_length']).mean()
print(df_qlen)
