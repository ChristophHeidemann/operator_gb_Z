# coding: utf-8
r"""
Ambiguities as defined for noncommutative Gröbner basis computations

Ambiguities characterise situations where one term can be reduced in
two different ways. The goal of noncommutative Gröbner basis computations
is to resolve all ambiguities by reducing S-polynomials formed from them. 

AUTHORS:

- Clemens Hofstadler (2023-03-01): initial version

"""

#############################################################################
#  Copyright (C) 2020 Clemens Hofstadler (clemens.hofstadler@jku.at).       #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 2, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from .auxiliary import flatten
from itertools import combinations, product
import os

############################################################################
# ambiguities
############################################################################
class Ambiguity:
    def __init__(self,ABC,Ai,Ci,Aj,Cj,i,j):
        self.__ABC = ABC
        self.__Ai = Ai
        self.__Ci = Ci
        self.__Aj = Aj
        self.__Cj = Cj
        self.__i = i
        self.__j = j
############################################################################
    def ABC(self): return self.__ABC
    def AC(self): return self.__ABC[:self.__Ai], self.__ABC[self.__Ci:],self.__ABC[:self.__Aj], self.__ABC[self.__Cj:]
    def ij(self): return self.__i, self.__j
    def i(self): return self.__i
    def j(self): return self.__j
    def degree(self): return len(self.__ABC)
############################################################################     
    def __eq__(self,other):
       return self.__i == other.__i and \
              self.__j == other.__j and \
              self.__Ai == other.__Ai and \
              self.__Aj == other.__Aj and \
              self.__Ci == other.__Ci and \
              self.__Cj == other.__Cj and \
              self.__ABC == other.__ABC
 ############################################################################                  
    def __ne__(self,other): return not (self == other)
#############################################################################
    def __hash__(self):
       return hash( (self.__ABC,self.__Ai,self.__Ci,self.__Aj,self.__Cj,self.__i,self.__j) )
#############################################################################
    def __repr__(self):
            Ai,Ci,Aj,Cj = self.AC()
            return "(" + self.__ABC + ", " + Ai + ", " + Ci + ", " + Aj + ", " + Cj + ", (" + str(self.__i) + ", " + str(self.__j) + "))"
############################################################################
    def __truediv__(self,other):
        r"""
        Divide ``self`` by ``other`` if possible.
        
        Am ambiguity `(ABC, A_i, C_i, A_j, C_j, i ,j)` is divisible by 
        another ambiguity `(A'B'C', A_i', C_i', A_j', C_j', i', j')` if
        there exist `L` and `R` such that `A_j = L A_j' and `C_j = C_j' R`.
        
        INPUT:
        
        - ``other`` -- Ambiguity
        
        OUTPUT:
        
        - `0` if ``other`` does not divides ``self``
        
        - `1` if ``other`` divides ``self`` and `L = R = 1`
        
        - `2` if if ``other`` divides ``self`` and `LR \neq 1`, i.e, ``other``
        properly divides ``self`` 
        
        TESTS::
        
            sage: a = Ambiguity()
        
        
        """  
        
        sAj, sCj = self.__Aj, self.__Cj
        oAj, oCj = other.__Aj, other.__Cj
        sdeg = self.degree()
        odeg = other.degree()
        
        start = sAj - oAj
        end = start + odeg
        
        
        
        if start < 0 or end > sdeg: return 0
        elif self.__ABC[start:end] == other.__ABC:
            if sdeg == odeg: return 1
            else: return 2
############################################################################    
    @staticmethod
    def generate_incls(a, b):        
        amb = []
        
        i,v = a
        j,w = b
        
        if i == j: return amb
           
        k = v.find(w, 0)
        while k >= 0:
            amb.append( Ambiguity(v,0,len(v),k,k+len(w),i,j) )
            k = v.find(w,k+1)

        return amb
            
############################################################################
    @staticmethod
    def generate_with_tries(prefix_trie, suffix_trie, words, i, amb_old, criterion=False):
        amb = []
        m = words[i]
        m_rev = m[::-1]
        len_m = len(m)
        for k in range(1,len_m):
            # overlap ambiguities with m = AB
            A = m[:-k]
            B = m[-k:]
            #amb += [Ambiguity(A+BC,len(A),len(A+BC),0,len_m,j,i) for j,BC in prefix_trie.values(B) if len(BC) > k and j <= i]
            for B_list, some_data in prefix_trie.values(B):
                for j,BC in B_list:
                    if len(BC) > k and j <=i:
                        amb += [Ambiguity(A+BC,len(A),len(A+BC),0,len_m,j,i)]

            # overlap ambiguities with m = BC
            B = m_rev[-k:]
            C = m[k:]            
            #amb += [Ambiguity(AB+C,0,len(AB),len(AB)-k,len(AB+C),j,i) for j,AB in suffix_trie.values(B) if len(AB) > k and j < i]
            for B_list, some_data in suffix_trie.values(B):
                for j,AB in B_list:
                    if len(AB) > k and j < i:
                        amb += [Ambiguity(AB+C,0,len(AB),len(AB)-k,len(AB+C),j,i)]


        # inclusion ambiguities with m = ABC       
        #amb += [Ambiguity(m,k-len(B)+1,k+1,0,len(m),j,i) for k,(j,B) in prefix_trie.iter(m) if j < i]
        for k,(B_list,some_data) in prefix_trie.iter(m):
            for j,B in B_list:
                if j < i:
                    amb += [Ambiguity(m,k-len(B)+1,k+1,0,len(m),j,i)]


        # inclusion ambiguities with m = B
        v = i,m
        #amb += flatten([Ambiguity.generate_incls((j,w),v) for j,w in enumerate(words[:i]) if len(w) > len_m])
        incl_ambs = []
        for j,w in enumerate(words[:i]):
            if len(w) < len_m:
                incl_ambs += [Ambiguity.generate_incls((j,w),v)]

        amb += flatten(incl_ambs)
                
        return amb
    

############################################################################
    @staticmethod       
    def generate_external_ambs(words, parent, maxdeg, oldlen):
        r"""
        Generate all external ambiguities for every pair of elements in G up to a certain max degree
        
        ........................
        
        INPUT:
        
        - ``G`` -- set of polynomials

        - ``parent`` -- free algebra where elements of G live in

        - ``maxdeg`` -- maximum length of leading monomials of external ambiguities

        - ``oldlen`` -- number of GB elements for which ambs were generated before
        
        OUTPUT:
        
        - `ex_ambs` -- external ambiguities up to certain degree 
        
        TESTS::
        
            sage: a = Ambiguity()
        
        
        """ 
        #check wether there are at least 2 elements in words 
        if len(words) < 2: return []

        vars = parent.translator().internal_gens() #translated variables of the free algebra
        pairs = list(combinations(words, 2)) #pairs of words


        #determine the max degree external overlaps can have
        sorted_pairs = sorted(pairs, key=len)
        max_ex_ov_deg = maxdeg - len(sorted(words, key=len)[0]) - len(sorted(words, key=len)[1])
        if max_ex_ov_deg < 0: return []

        #generate external overlaps
        ex_overlaps = ['']
        for i in range(max_ex_ov_deg):
            combs = product(vars, repeat=i+1)
            ex_overlaps.extend([''.join(c) for c in combs])

        #generate external ambiguities by iterating over pairs of words
        #we are interestet in the ex amb ABC and CBA
        ex_ambs = []

        #this loop is for the ex ambs between the new GB elements
        for i in range(oldlen, len(words)):
            for j in range(i+1, len(words)):
                A = words[i]
                C = words[j]
                for B in ex_overlaps:        
                    if len(A + B + C) > maxdeg: break
                    ex_ambs += [Ambiguity(A+B+C, 0, len(C), len(A) + len(B), len(A+B+C), i, j)]
                    ex_ambs += [Ambiguity(C+B+A, 0, len(A), len(C) + len(B), len(C+B+A), j, i)] 

        #this loop is for the ex ambs between a new and an old GB element
        for i in range(oldlen):
            for j in range(oldlen, len(words)):
                A = words[i]
                C = words[j]
                for B in ex_overlaps:        
                    if len(A + B + C) > maxdeg: break
                    ex_ambs += [Ambiguity(A+B+C, 0, len(C), len(A) + len(B), len(A+B+C), i, j)]
                    ex_ambs += [Ambiguity(C+B+A, 0, len(A), len(C) + len(B), len(C+B+A), j, i)] 

        #amb += [Ambiguity(AB+C,0,len(AB),len(AB)-k,len(AB+C),j,i) for j,AB in suffix_trie.values(B) if len(AB) > k and j < i]
        return ex_ambs

############################################################################       
    def shrink(self):
    
        ABC = self.ABC()
        i,j = self.ij()
        Ai,Ci,Aj,Cj = self.__Ai, self.__Ci, self.__Aj, self.__Cj
        k = min(Ai,Aj)
        l = max(Ci,Cj)
        
        return Ambiguity(ABC[k:l], Ai-k, Ci-k, Aj-k, Cj-k,i,j) 
