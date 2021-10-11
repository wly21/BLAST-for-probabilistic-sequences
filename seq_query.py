import random
import os

def read_sequence(demo_seq_len = 0):
    '''

    Parameters
    ----------
    demo_seq_len : int, optional
        Since the entire Boreo Eutherian Sequence is long (length = 604466), 
        we provide an option to use only part of the sequence for testing and demonstration.
        The first demo_seq_len characters (nucleotides), with corresponding data, will be returned

    Returns
    -------
    SEQUENCE : list of characters, the most probable nucleotide at each position
    SEQ_PROB : list of float numbers, the probability of the the most probable nucleotide at each position
    NEUC_PROB: list of dict, the probabilities for different nucleotides at each position
    Q_PROB: dict, the overall probability of each nucleotide in the entire sequence, used in mutate() 

    '''    
    
    dirname = os.path.dirname(__file__)
    seq_file = os.path.join(dirname, 'data/BoreoEutherian_sequence.txt')
    prob_file = os.path.join(dirname, 'data/BoreoEutherian_prob.txt')
    
    SEQUENCE = []
    with open(seq_file,'r') as f:
        for line in f.readlines():
            for neuc in line.strip():
                if neuc!=' ':
                    SEQUENCE.append(neuc)
    
    SEQ_PROB = []
    with open(prob_file,'r') as f:
        for line in f.readlines():
            SEQ_PROB = line.strip().split(' ')
            SEQ_PROB = [ float(prob) for prob in SEQ_PROB]
            
    if(demo_seq_len != 0):
        SEQUENCE = SEQUENCE[:demo_seq_len]
        SEQ_PROB = SEQ_PROB[:demo_seq_len]
    
    NEUC_PROB = [{'A':0,'C':0,'G':0,'T':0} for i in range(len(SEQ_PROB))]
    NEUC = ['A','C','G','T']
    for i in range(len(SEQ_PROB)):
        NEUC_PROB[i][SEQUENCE[i]] = SEQ_PROB[i]
        for neucl in NEUC:
            if NEUC_PROB[i][neucl]==0:
                NEUC_PROB[i][neucl] = (1-SEQ_PROB[i])/3
                
    Q_PROB = {'A':0,'C':0,'G':0,'T':0}
    for i in range(len(SEQ_PROB)):
        for neucl in NEUC:
            Q_PROB[neucl] += NEUC_PROB[i][neucl]
    for neucleotide in NEUC:
        Q_PROB[neucleotide] /= len(SEQ_PROB)
        
    return SEQUENCE,SEQ_PROB, NEUC_PROB, Q_PROB


def generate_query(SEQUENCE,SEQ_PROB,query_len, sub_rate=0.1, del_rate=0.02, ins_rate=0.02):
    # generate query based on the probability provided
    start_pos = random.randint(0, len(SEQUENCE) - query_len)
    query = ""
    for i in range(query_len):
        p = random.random()
        if p <= SEQ_PROB[start_pos+i]:
            query += SEQUENCE[start_pos+i]
        else:
            possible_nucleotide = ['A','C','G','T']
            possible_nucleotide.remove(SEQUENCE[start_pos+i])
            query += possible_nucleotide[random.randint(0,2)]            
    query = mutate(query, sub_rate, del_rate, ins_rate)
    return start_pos, query
    
def mutate(query, sub_rate, del_rate, ins_rate):
    # mutate based on error rate
    result = ""
    for s in query:
        if random.random() <= del_rate:
            #delete
            continue
        else:
            if random.random() <= sub_rate:
                #substitute
                possible_nucleotide = ['A','C','G','T']
                possible_nucleotide.remove(s)
                result += possible_nucleotide[random.randint(0,2)]
            else:
                #no change
                result += s           
        if random.random() <= ins_rate:
            #insert
            possible_nucleotide = ['A','C','G','T']
            result += possible_nucleotide[random.randint(0,3)]
            if random.random() <= ins_rate*ins_rate:
                #second insertion
                possible_nucleotide = ['A','C','G','T']
                result += possible_nucleotide[random.randint(0,3)]
    return result