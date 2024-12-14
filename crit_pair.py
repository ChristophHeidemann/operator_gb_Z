from __future__ import absolute_import

#from .cofactor import Cofactor

from copy import copy
from sage.all import xgcd # type: ignore
from sage.all import lcm # type: ignore

############################################################################
# S-polynomials
############################################################################
class CritPair:
    def __init__(self, amb , gi, gj):
    
        A = gi.parent()
        
        self.__degree = amb.degree()
        Ai,Ci,Aj,Cj = amb.AC()
        g1 = gi.lrmul(Ai,Ci)
        g2 = gj.lrmul(Aj,Cj)
        #Bezout coefficients for G-Polynomials
        u = xgcd(gi.lc(),gj.lc())[1]
        v = xgcd(gi.lc(),gj.lc())[2]
        #lcm multipliers for S-Polynomials
        qf = lcm(gi.lc(),gj.lc()) / gi.lc()
        qg = lcm(gi.lc(),gj.lc()) / gj.lc() 



        if g1 == g2:
            self.__degree = -1   
            self.__f = g1
            self.__g = g2
            #Für G-Pols
            self.__u = u
            self.__v = v

        else:
            #i,j = amb.ij()
            #g1.append_cofactor(Cofactor(1,Ai,i,Ci,A))  #Chris: 'cofactor ausgeklammert'
            #g2.append_cofactor(Cofactor(1,Aj,j,Cj,A))
            
            #Für S-Pols
            self.__f = g1
            self.__g = g2
            self.__qf = qf
            self.__qg = qg

            #Für G-Pols
            self.__u = u
            self.__v = v
            
############################################################################
    def degree(self): return self.__degree
    def f(self): return self.__f
    def g(self): return self.__g
    def fg(self): return [self.__f.coeff_mul(self.__qf), self.__g.coeff_mul(self.__qg)]
    def uv(self): return [self.__f.coeff_mul(self.__u), self.__g.coeff_mul(self.__v)] #Critical pair for G-Pol
############################################################################
    def __repr__(self):
        return "(" + str(self.__f) + ", " + str(self.__g) + ")"
