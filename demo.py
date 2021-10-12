import BLAST_prob
from BLAST_prob import seq_query, BLAST

# read sequence data
SEQUENCE,SEQ_PROB, NEUC_PROB, Q_PROB = seq_query.read_sequence(demo_seq_len=2500)

# generate query for test/demo
start_pos, query = seq_query.generate_query(SEQUENCE,SEQ_PROB,query_len=100)

# query information
print("The query is taken from sequence position " + str(start_pos) + " to " + str(start_pos+100-1))

# align the query to the probablistic sequence
result_hit, t = BLAST.BLAST(NEUC_PROB,query)
print("result alignment:", result_hit)
print("time by step:",t)
