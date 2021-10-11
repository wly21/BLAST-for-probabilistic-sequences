from BLAST_prob import Trie, Hit

def BLAST(nucl_prob, query, word_l = 9, score_thres = 6, window = 30):
    print("Entered BLAST////////////////////////////////////////")
    TRIE = Trie.possible_matching(query,word_l, score_thres)
    print("Trie created")
    ONEHIT = Hit.find_one_hit(TRIE, nucl_prob, word_l, query,score_thres/2)
    print("One hit found")
    TWOHIT = Hit.find_two_hit(nucl_prob,query,ONEHIT,window,score_thres/2)
    print("two hit found")
    UNGAP_HIT = Hit.ungapped_Extension(TWOHIT, nucl_prob, query)
    print("hsp found")
    GAPPED_HIT = Hit.gapped_extension(UNGAP_HIT, nucl_prob, query)
    print("gapped extension done")
    MAX_HIT = Hit.find_max_hit(GAPPED_HIT)
    print("Max hit found")
    
    return MAX_HIT



