# coding: utf-8
"""
Implementation of noncommutative ideals in free algebras

AUTHOR:

- Clemens Hofstadler (2023-03-01): initial version

"""

# ****************************************************************************
#                          Copyright (C) 2023
#      Clemens Hofstadler(clemens.hofstadler@mathematik.uni-kassel.de)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import itertools
from copy import deepcopy

from .f4 import F4
from .free_algebra import MyFreeAlgebra
from .nc_polynomial import NCPolynomial
#from .nc_ideal_right import NCIdeal_right
#from .normal_form import reduced_form, interreduce
#from .quiver import Quiver

############################################################################
#  NCIdeal
############################################################################
class NCIdeal:
    def __init__(self,gens,order=None):
        
        F = gens[0].parent()
        X = order if order else F.gens()
        A = MyFreeAlgebra(F.base_ring(),X)        
        
        self.__gens = gens
        self.__internal_gens = [f if isinstance(f,NCPolynomial) else A(f) for f in gens]
        self.__parent = A
        self.__G = []
        self.__algo = None
############################################################################
    def gens(self): return self.__gens
    def parent(self): return self.__parent
    def internal_gens(self): return self.__internal_gens
    def order(self): return self.__parent.order()
    def algo(self): return self.__algo
    def base_ring(self): return self.__parent().base_ring()
############################################################################    
    def __repr__(self):
        return "NCIdeal %s of %s" % (str(tuple(self.__gens)),str(self.__parent))
############################################################################
    def groebner_basis(self,maxdeg,maxiter=10,trace_cofactors=True,criterion=True,reset=True,verbose=0):
        if not self.__algo:
            self.__algo = F4(self)
        self.__G = self.__algo.compute_basis(maxiter=maxiter,maxdeg=maxdeg,trace_cofactors=trace_cofactors,criterion=criterion,reset=reset,verbose=verbose)
        return self.__G
############################################################################
    def __add__(self,other):
        if self.__parent != other.__parent:
            raise ValueError("Ideals have to be defined over the same ring.")
        
        return NCIdeal(self.__gens + other.__gens, order=self.__order) 
############################################################################






