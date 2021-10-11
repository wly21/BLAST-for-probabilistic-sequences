import itertools
#import torch
from BLAST_prob import Trie

class Hit:
    def __init__(self, dstart=None,dend=None, qstart=None,qend=None,score=None,diag=None):
        self.db_start = dstart
        self.db_end = dend
        self.query_start = qstart
        self.query_end = qend
        self.prob_score = score
        self.diagonal = diag
    
    def __str__(self):
        return "The hit matches "+str(self.db_start)+" and "+str(self.db_end)+" in the database sequence with the positions from "+str(self.query_start)+" to "+str(self.query_end)+" in query sequence with score "+str(self.prob_score)
    

def find_one_hit(search_trie, sequence_prob, word_len, query, score_threshold):
    NEUC = ['A','C','G','T']
    seq_words = []
    hits = []
    
    for i in range(len(sequence_prob)-word_len+1):
        prob_words = []
        
        for t in range(word_len):
            prob_neucleotide = {k:v for (k,v) in sequence_prob[i+t].items() if v!=0}
            prob_words = [k for (k,v) in prob_neucleotide.items()]
            
            if t==0:
                seq_words = [[k] for (k,v) in prob_neucleotide.items()]
            else:
                seq_words = list(itertools.product(seq_words, prob_words))
                temp_seq = []
                for p in range(len(seq_words)):
                    tempword = list(seq_words[p][0])
                    tempword.append(seq_words[p][1])
                    temp_seq.append(tempword)
                seq_words = temp_seq
        
        for j in range(len(seq_words)):
            match_query = Trie.SearchTrie.Detect_Hit(search_trie, seq_words[j])
            if match_query is None:
                continue
            else:
                for query_w in match_query:
                    qword = query[query_w[0]:query_w[0]+word_len]
                    score = 0
                    for p in range(word_len):
                        if qword[p]==seq_words[j][p]:
                            score += sequence_prob[i+p][seq_words[j][p]] * 1
                        else:
                            score -= sequence_prob[i+p][seq_words[j][p]] * 2
                    if score>=score_threshold:
                        hits.append(Hit(i,i+word_len-1,query_w[0],query_w[0]+word_len-1,score,None))
        if i>0 and i%100==0:
            print("One-hit: Checking {} positions done".format(i))                 
    return hits
       
#two-hit, joining one-hits
def find_two_hit(sequence_prob,query,hits,window,score_threshold):
    two_hits=[]
    for i in range (0,len(hits)):
        hits[i].diagonal = hits[i].db_start - hits[i].query_start  #db_pos-q_pos
    for i in range (0,len(hits)-1):
        hit1 = hits[i]
        j=i+1
        hit2 = hits[j]
        while hit2.db_end < (hit1.db_start + window):
            if(hit1.diagonal==hit2.diagonal and hit1.db_end<hit2.db_start):
                #calculate total score
                inter_db = sequence_prob[hit1.db_end+1 : hit2.db_start]
                inter_q = query[hit1.query_end+1 : hit2.query_start]
                
                inter_score = 0
                for k in range(0,len(inter_db)):
                    s=1*inter_db[k][inter_q[k]] + (-2)*(1-inter_db[k][inter_q[k]])
                    inter_score += s
                
                score = hit1.prob_score + hit2.prob_score + inter_score
                if(score>score_threshold):
                    #create new hit and append to two_hits
                    h = Hit(hit1.db_start,hit2.db_end,hit1.query_start,hit2.query_end,score)
                    two_hits.append(h)
            j+=1
            if j<len(hits):
                hit2 = hits[j]
            else:
                break
    return two_hits 


def ungapped_extension(hits, sequence_prob, query):
    ungapped_hits = []
    for hit in hits:
        ext_hit = Hit(hit.db_start,hit.db_end,hit.query_start,hit.query_end,hit.prob_score,hit.diagonal)
        
        # first extend on the right of the second hit
        acc_score = 0
        db_e = hit.db_end
        qy_e = hit.query_end
        for i in range(qy_e+1,len(query)):
            if db_e+(i-qy_e)>=len(sequence_prob):
                ext_hit.db_end = db_e+(i-qy_e)-1
                ext_hit.query_end = i-1
                break
            prob_neucs = {k:v for (k,v) in sequence_prob[db_e+(i-qy_e)].items() if v!=0}
            prob_words = [k for (k,v) in prob_neucs.items()]
            tempscore = acc_score
            if query[i] in prob_words:
                acc_score += prob_neucs[query[i]]*1
                prob_words.remove(query[i])
            for k in prob_words:
                acc_score -= prob_neucs[k]*2
            # stop if the accumulated score becomes lower than 0
            if acc_score<0:
                ext_hit.db_end = db_e+(i-qy_e)-1
                ext_hit.query_end = i-1
                acc_score = tempscore
                break
            if i==len(query)-1:
                ext_hit.db_end = db_e+(i-qy_e)
                ext_hit.query_end = i
        ext_hit.prob_score = hit.prob_score + acc_score
        
        # then extend on the left of the first hit
        acc_score = 0
        db_s = hit.db_start
        qy_s = hit.query_start
        for i in range(qy_s-1,-1, -1):
            if db_s-(qy_s-i)<0:
                ext_hit.db_start = db_s-(qy_s-i)+1
                ext_hit.query_start = i+1
                break
            prob_neucs = {k:v for (k,v) in sequence_prob[db_s-(qy_s-i)].items() if v!=0}
            prob_words = [k for (k,v) in prob_neucs.items()]
            tempscore = acc_score
            if query[i] in prob_words:
                acc_score += prob_neucs[query[i]]*1
                prob_words.remove(query[i])
            for k in prob_words:
                acc_score -= prob_neucs[k]*2
            # stop if the accumulated score becomes lower than 0
            if acc_score<0:
                ext_hit.db_start = db_s-(qy_s-i)+1
                ext_hit.query_start = i+1
                acc_score = tempscore
                break
            if i==0:
                ext_hit.db_start = db_s-(qy_s-i)
                ext_hit.query_start = i
        ext_hit.prob_score += acc_score
        
        ungapped_hits.append(ext_hit)

    return ungapped_hits

# Gapped extension: Smith-Waterman algorithm modified for the probabilistic sequence
def gapped_extension(hits, sequence_prob, query):
    gapped_hits = []
    #dtype = torch.cuda.FloatTensor 
    for hit in hits:
        ext_hit = Hit(hit.db_start,hit.db_end,hit.query_start,hit.query_end,hit.prob_score,hit.diagonal)
        
        if hit.query_end!=len(query)-1:
            M_r = [[0 for i in range(len(query)-hit.query_end)] for j in range(len(query)-hit.query_end)]
            STATES = ['MATCH','GAP_IN_DB','GAP_IN_QUERY']
            best_score = 0
            for i in range(len(query)-hit.query_end):
                for j in range(len(query)-hit.query_end):
                    if i==0 and j==0:
                        continue
                    if i == 0 and j>0:
                        M_r[i][j] = M_r[i][j-1] -2
                    elif j == 0 and i>0:
                        M_r[i][j] = M_r[i-1][j] -2
                    else:
                        if hit.db_end+j<len(sequence_prob):
                            prob_neucs = {k:v for (k,v) in sequence_prob[hit.db_end+j].items() if v!=0}
                        else:
                            prob_neucs = {}
                        match_score = M_r[i-1][j-1]
                        for k,v in prob_neucs.items():
                            if k==query[hit.query_end+i]:
                                match_score += v*1
                            else:
                                match_score -= v*2
                        qy_gap_score = M_r[i-1][j] -2
                        db_gap_score = M_r[i][j-1] -2

                        if match_score >= qy_gap_score and match_score >= db_gap_score:
                            M_r[i][j] = match_score
                        elif db_gap_score > match_score and db_gap_score > qy_gap_score:
                            M_r[i][j] = db_gap_score
                        else:
                            M_r[i][j] = qy_gap_score

            best_score = max(M_r[-1][:])
            best_j_index = M_r[-1][:].index(best_score)
            ext_hit.query_end = len(query)-1
            ext_hit.db_end = min(ext_hit.db_end+best_j_index, len(sequence_prob)-1)
            ext_hit.prob_score += best_score

        best_score = 0
        if hit.query_start!=0:
            M_l = [[0 for i in range(hit.query_start)] for j in range(hit.query_start)]
            for i in range(hit.query_start):
                for j in range(hit.query_start):
                    if i==0 and j==0:
                        continue
                    if i == 0 and j>0:
                        M_l[i][j] = M_l[i][j-1] -2
                    elif j == 0 and i>0:
                        M_l[i][j] = M_l[i-1][j] -2

                    else:
                        if hit.db_end+j<len(sequence_prob):
                            prob_neucs = {k:v for (k,v) in sequence_prob[hit.db_start-j].items() if v!=0}
                        else:
                            prob_neucs = {}
                        match_score = M_l[i-1][j-1]
                        for k,v in prob_neucs.items():
                            if k==query[hit.query_start-i]:
                                match_score += v*1
                            else:
                                match_score -= v*2
                        qy_gap_score = M_l[i-1][j] -2
                        db_gap_score = M_l[i][j-1] -2

                        if match_score >= qy_gap_score and match_score >= db_gap_score:
                            M_l[i][j] = match_score
                        elif db_gap_score > match_score and db_gap_score > qy_gap_score:
                            M_l[i][j] = db_gap_score
                        else:
                            M_l[i][j] = qy_gap_score

            best_score = max(M_l[-1][:])
            best_j_index = M_l[-1][:].index(best_score)
            ext_hit.query_start = 0
            ext_hit.db_start = max(ext_hit.db_start-best_j_index, 0)
            ext_hit.prob_score += best_score
        gapped_hits.append(ext_hit)
    return gapped_hits


def find_max_hit(hits):
    max_hit = None
    max_score = None
    for hit in hits:
        if max_score is None or hit.prob_score>max_score:
            max_hit = hit
            max_score = hit.prob_score
    
    return max_hit