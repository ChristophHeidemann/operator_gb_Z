from time import time
from copy import copy

import ahocorasick

from .ambiguity import *
#from .cofactor import Cofactor
#from .free_algebra import MyFreeAlgebra
from .nc_monomial import NCMonomial
from .nc_polynomial import NCPolynomial
#from .normal_form import interreduce
from .crit_pair import CritPair
#from .linear_algebra import faugere_lachartre
from sage.all import matrix, ZZ #type: ignore

class F4():

    def __init__(self,I):
        self.__parent = I.parent()
        self.__gens = I.gens()
        self.__internal_gens = I.internal_gens()
        self.__criterion = True
        self.__maxdeg = -1
        self.__lm = ahocorasick.Automaton()
        self.__suffix_trie = ahocorasick.Automaton()
        self.__verbose = 0
        self.__amb = []
        self.__G = []
        self.__constant_flag = False
############################################################################        
    def lm(self): return list(self.__lm.items()) 
    def parent(self): return self.__parent
    def gens(self): return self.__gens
    def internal_gens(self): return self.__internal_gens
    def amb(self): return self.__amb
    def G(self): return self.__G
############################################################################
    def clear(self):
        self.__lm.clear()
        self.__suffix_trie.clear()
        self.__amb = []
        self.__G = []
        self.__constant_flag = False
############################################################################
# Compute Gröbner basis
############################################################################
    def compute_basis(self,maxiter=10,maxdeg=-1,trace_cofactors=True,criterion=False,reset=True,verbose=0):
        global nr_pairs
        global zero_reductions
        
        self.__maxdeg = maxdeg
        self.__criterion = criterion
        self.__verbose = verbose
        
        nr_pairs = 0
        zero_reductions = 0
        
        if reset or not self.__G: 
            self.prepare_input(verbose,trace_cofactors=trace_cofactors)
            # constant_flag checks if GB contains 1
            if self.__constant_flag: return self.__G
            self.compute_ambiguities()
        

        self.__already_rewritten = len(self.__G)
        oldlen = len(self.__G)    
        count = 0
        maxiter = 10

        while count < maxiter and self.__amb:
            count += 1

            PP = self.reduction()

            #check whether new leading terms came from the reduction and add them
            PP = self.check_leading_terms(PP)

            self.add_polynomials(PP)

            if self.__constant_flag: 
                print('GB contains 1')
                break            

            self.compute_ambiguities(oldlen)
            oldlen = len(self.__G)
            


        print('Gröbner basis: ' + str(self.__G))
        return self.__G

############################################################################
    def compute_ambiguities(self,oldlen=0):
        
        G = self.__G
        parent = self.__parent
        maxdeg = self.__maxdeg
        verbose = self.__verbose
        criterion = self.__criterion
                
        words = [str(g.lm()) for i,g in enumerate(G)]
        l = len(self.__amb)
            
        # compute ambiguities    
        start = time()
        prefix_trie = self.__lm
        suffix_trie = self.__suffix_trie
        for i in range(oldlen,len(words)):
            amb_i = Ambiguity.generate_with_tries(prefix_trie,suffix_trie,words,i,self.__amb,criterion=criterion)
            if maxdeg > 0: amb_i = [a for a in amb_i if a.degree() <= maxdeg]
            #if criterion:
                #self.__amb = Ambiguity.chain_criterion(self.__amb,amb_i,words[i],i)
            self.__amb += amb_i
                                    
        #generate external ambiguities
        if oldlen != len(words):
            self.__amb += Ambiguity.generate_external_ambs(words, parent, maxdeg, oldlen)


        if verbose > 0:
            print("%d ambiguities in total (computation took %.5f)" % (len(self.__amb), time()-start))
                    
        self.__amb.sort(key = lambda p : p.degree())


############################################################################
# Reduction & Symbolic Preprocessing
############################################################################
    def reduction(self,trace_cofactors=False):  #Chris: 'auf False gesetzt'
        global nr_pairs
        global zero_reductions
                
        G = self.__G
        amb = self.__amb
        verbose = self.__verbose
        
        # choose ambiguities            
        d = amb[0].degree()
        P = [a for a in amb if a.degree() == d] 
        self.__amb = amb[len(P):]       #remove ambiguities that are processed
        P = [CritPair(a,G[a.i()],G[a.j()]) for a in P]        #Generate parts for crit pairs
        #copy list so G-pols are complete
        P_G = P
        # remove zero S-polynomials
        P = [pair for pair in P if pair.degree() > -1]
        if not P: return []
        if verbose > 0:
            print("%d critical pairs will be reduced." % len(P))
        
        nr_pairs += len(P) 

        # do symbolic preprocessing
        F = [pair.fg() for pair in P]
        F = flatten(F)                  #Critical pairs for S-Pols
        F_G_pol = [pair.uv() for pair in P_G]
        F_G_pol = flatten(F_G_pol)            #Critical pairs for G-Pols
        F += F_G_pol                                             #joined list of crit pairs
        F = list(set(F))                    #remove duplicate elements
        F = [f for f in F if f.lc()!=0]         #remove zeros        
    
        pivot_rows,pivot_columns,columns = self.symbolic_preprocessing(F)

        columns = list(columns)
        columns.sort(reverse=True)              #sort the appearing monomials 
        columns = {m:i for (i,m) in enumerate(columns)}

        F += pivot_rows             #assemble the matrix M_{P U G'}
        M = self.set_up_matrix(F,columns)               #set up the Matrix M_{P U G'}
        M = M.hermite_form()

        PP = self.matrix_to_polies(M,columns)             #translate the matrix back into polynomial form

        zero_reductions += len(P) - len(PP) 
            
        return PP
############################################################################
    def symbolic_preprocessing(self,F):
        done = set()
        todo = {m for f in F for m in f.monomials()}

        reducers = []
        pivot_columns = set()
        while todo:
            m = todo.pop()
            done.add(m)
            agb = self.find_reducer(str(m))
            if agb:
                assert m == agb.lm()
                reducers.append(agb)
                pivot_columns.add(agb.lm())
                todo.update(set(agb.monomials()).difference(done))
        return reducers,pivot_columns,done
############################################################################
    def find_reducer(self,m):
        lm = self.__lm
        G = self.__G
                
        reducer = list(lm.iter(m))
        if not reducer: return None
        k,(l,data) = min(reducer, key=lambda r:r[1][0]) # choose red (degree)
        (i,lm) = l[0]
        g = G[i]
        a = m[:k-len(lm)+1]
        b = m[k+1:]
        agb = g.lrmul(a,b)
        return agb

############################################################################        
# Auxiliary
############################################################################      
    def prepare_input(self,verbose,trace_cofactors=True):
        
        self.clear()
        # change data structure
        A = self.__parent
        G = [copy(f) for f in self.__internal_gens]
        
        # make monic and add cofactors
        #for i,g in enumerate(G):
            #c = g.make_monic()
            #g.reset_cofactors()
            #if trace_cofactors: g.append_cofactor(Cofactor(1/c,'',i,'',A))
        
        # interreduce generators
        #G = interreduce(G)
        #if verbose > 0:
            #print("Interreduced the generators from %d elements to %d elements.\n" % (len(self.__gens),len(G)))
        
        # add interreduced generators to GB    
        self.add_polynomials(G)

############################################################################                    
    def add_polynomials(self,PP):
        oldlen = len(self.__G)
        for i,p in enumerate(PP):
            m = str(p.lm())
            if m:
                self.add_new_key(self.__lm, m, (i+oldlen,m), 'some data')
                #self.add_new_key(self.__lm, m, i+oldlen, m)
                
                #self.__lm.add_word(m,(i+oldlen,m))          #ersetzt durch die obere zeile
                #self.__suffix_trie.add_word(m[::-1],(oldlen+i,m))
                self.add_new_key(self.__suffix_trie, m[::-1],(oldlen+i,m), 'some data')
            else:
                self.__constant_flag = True
                break

        self.__lm.make_automaton()
        self.__G += PP

############################################################################        
    @staticmethod
    def add_new_key(A,new_key,new_idx,new_data):
        idxs, data = [], new_data
        if A.exists(new_key):
            idxs, data = A.get(new_key)
        A.add_word(new_key,(idxs + [new_idx], data))
############################################################################                    
    @staticmethod
    def set_up_matrix(rows, columns): 
        # Determine the number of rows and columns
        nr = len(rows)
        nc = len(columns)
        
        # Create a sparse integer matrix with the given number of rows and columns
        A = matrix(ZZ, nr, nc, sparse=True)
        
        # Fill the matrix with coefficients from rows
        for i, f in enumerate(rows):
            for c, m in zip(f.coefficients(), f.monomials()):
                j = columns[m]
                A[i, j] = c  # Set the value in the matrix A at position (i, j)
        
        return A
############################################################################        
    def matrix_to_polies(self,M,columns):
        columns = list(columns)
        #compute polynomials
        pos = M.nonzero_positions()
        if not pos: return []
        rank = pos[-1][0]+1
        coefficients = [[] for i in range(rank)]
        monomials = [[] for i in range(rank)]
        for (i,j),c in reversed(list(M.dict().items())):            
            coefficients[i].append(c)
            monomials[i].append(columns[j])

        PP = [NCPolynomial(c,m) for c,m in zip(coefficients,monomials)]

        PP.sort(key=lambda f : f.lm())
        return PP   
############################################################################        
    def check_leading_terms(self, PP):
        new_PP = []
        for p in PP:
            contained = False
            for g in self.__G:
                if p.lm() == g.lm():
                    if p.lc()%g.lc() == 0:
                        contained = True
                        break
            if not contained: new_PP.append(p) 

        return new_PP