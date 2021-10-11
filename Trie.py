from functools import reduce
import operator

def Score_Between_Words(word1, word2):
    # We applied the scoring scheme (Match: +1, Mismatch: -2, Gap Costs: linear) based on the suggestion on the BLAST website. 
    # Read more here: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp#Reward-penalty
    score = 0
    for i in range(len(word1)):
        if word1[i]==word2[i]:
            score += 1
        else:
            score -= 2
    
    return score

# generate a search trie such that each leaf represents a possible word that has a match in query
def possible_matching(query, word_len, threshold):
    candidates = [['A'],['C'],['G'],['T']]
    # generate a list of all possible words (i.e. letter/nucleotide permutations) with given word length
    for i in range(1,word_len):
        temp_cand = []
        for cand in candidates:
            cand.append('A')
            temp_c = list(cand[:-1])
            temp_c.append('C')
            temp_cand.append(temp_c)
            temp_c = list(cand[:-1])
            temp_c.append('G')
            temp_cand.append(temp_c)
            temp_c = list(cand[:-1])
            temp_c.append('T')
            temp_cand.append(temp_c)
        candidates.extend(temp_cand)
    print("Finding candidates done")
    
    # generate all words in query sequence with given word length
    query_words = []
    for j in range(len(query) - word_len +1):
        if query[j:j+word_len] not in query_words:
            query_words.append(query[j:j+word_len])
    print("Query words done")
    
    # filter the candidates with the given threshold
    qualified = [False for i in range(len(candidates))]
    for w in range(len(candidates)):
        for q in query_words:
            if Score_Between_Words(candidates[w],q)>threshold:
                qualified[w] = True
    
    stay_candidates = []
    for p in range(len(qualified)):
        if qualified[p]:
            stay_candidates.append(p)
    tempcand = [candidates[t] for t in stay_candidates]
    candidates = tempcand
    print("Qualified")
    
    print("Begin to build the search trie")
    
    # transfer the list of candidates into a search tree - a dictionary
    return SearchTrie(candidates, query_words, threshold, Score_Between_Words)
    
class SearchTrie:
    def __init__(self, candidates, query_words, threshold, Calculate_Score):
        self.trie = {}
        try:
          match = candidates[0]
        except IndexError:
          print("No possible matchings for hit")
        for fs in range(1,len(match)):
            temp_f = match[:fs]
            reduce(operator.getitem, temp_f[:-1], self.trie)[temp_f[-1]] = {}
        reduce(operator.getitem, match[:-1], self.trie)[match[-1]] = []

        print("There are",len(candidates),"possible matchings")
        for m in range(1,len(candidates)):
            match = candidates[m]
            prev = self.trie
            for fs in range(1,len(match)):
                temp_f = match[:fs]
                if temp_f[-1] not in prev.keys():
                    prev[temp_f[-1]] = {}
                prev = reduce(operator.getitem, temp_f[:-1], self.trie)[temp_f[-1]]
            reduce(operator.getitem, match[:-1], self.trie)[match[-1]] = []
        print("Branches created")

        for t in candidates:
            for index in range(len(query_words)):
                score = Calculate_Score(t,query_words[index])
                if score>threshold:
                    reduce(operator.getitem, t[:-1], self.trie)[t[-1]].append( (index, score) )

        #print(self.trie)
        print("Record scores done")
    
    def Detect_Hit(self, word):
        current_path = self.trie
        for i in range(len(word)):
            if word[i] not in current_path.keys():
                return
            tempword = word[:i+1]
            current_path = reduce(operator.getitem, tempword[:-1], self.trie)[tempword[-1]]
        return reduce(operator.getitem, word[:-1], self.trie)[word[-1]]
    
