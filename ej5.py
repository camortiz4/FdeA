# -*- coding: utf-8 -*-
"""
@author: Camilo Ortiz

Created on Wed Oct 21 10:05:11 2020

"""

import numpy as np
import math as math
import matplotlib
import matplotlib.pyplot as plt
from sympy.solvers import solve
from sympy import Symbol


sondeo = np.loadtxt('./sondeo.txt',skiprows=4,)

"""
PRES	HGHT	TEMP	DWPT	RELH	MIXR	DRCT	SKNT
hPa 	m   	C	    C	    %	    g/kg	deg 	knot
"""

l = 2.5e6
cp = 1004
rd = 287.05
rv = 461.5
k = rd/cp
T_tr = 273.15
e_tr = 611.73
epsilon = rd/rv


p_o= sondeo[0,0]*100
p = sondeo[:,0]*100
T = sondeo[:,2]+273
Td = sondeo[:,3]+273
w = sondeo[:,5]/1000




def potencial(T,P):
    return  T*(p_o/P)**k

def potencialEquivalente(T,P,w):
    return  T*(p_o/P)**(k)*math.exp(l*w/(cp*T))

def clausiusClapeyron(T):
    return math.exp(-l/rv*(T**-1 - T_tr**-1) + math.log(e_tr))

def mixingRatio(T,P):
    return epsilon*clausiusClapeyron(T)/P

def gammaDry(T_o,P_o,P):
    return potencial(T_o, P_o)*(P/p_o)**k

# def gammaMoist(T_par,P_par,w_par, p_o,T):
#     return T * p_o**k / potencialEquivalente(T_par, P_par, w_par) \
#            * math.exp(l*mixingRatio(T, P) /(cp * T) )

def gammaMoist(T_par,P_par,w_par, p_o,P):
    T = Symbol('T')
    return solve(T * (p_o/P)**k * math.exp(l*mixingRatio(T, P) /(cp * T))\
           - potencialEquivalente(T_par, P_par, w_par) 
           ,T)
    


tempPotencialEq = []

for i in range(len(p)):
    tempPotencialEq.append(potencialEquivalente(T[i], p[i],w[i]))
        

i_parcel = tempPotencialEq.index(max(tempPotencialEq[0:50]))
e_sat_parcela = clausiusClapeyron(Td[i_parcel])

NCA = epsilon*e_sat_parcela/w[i_parcel]


T_parcela = []
p_o = p[i_parcel]
for i in range(len(p)):
    T_parcela.append(gammaDry(T[i_parcel], p[i_parcel], p[i]))
        
print(gammaMoist(T[i_parcel], p[i_parcel], w[i_parcel], p_o, 50000))

fig, ax = plt.subplots()
ax.plot(T,p/100)
ax.plot(Td,p/100)
ax.plot(T_parcela ,p/100)
ax.set_xlim(min(T)-5,max(T)+5)
ax.set_ylim(1100,100)
plt.yscale('log')
ax.set_xlabel('T(K)')
ax.set_ylabel('log(P)')
ax.set_yticks(np.arange(100,1100,100))
ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

plt.show()



# e_parcela = []

# for i in range(len(p)):
#     e_parcela.append(w[index_parcela]*p[i]/epsilon)

# for i in range(len(e_parcela)-1):
#     if e_parcela[i] >= e_sat_parcela and e_sat_parcela >= e_parcela[i+1]:
#         NCA = p[i+1]/100
#         print(NCA)
#         index_NCA = i+1

# T_parcela = []

# for i in range(index_NCA):
#     T_parcela.append(potencial(T[index_parcela], p[index_parcela])\
#                      *(p[i]/p_o)**(k))
# for i in range(index_NCA,len(p)):
#     T_parcela.append(w[index_parcela]*l/cp * math.log(potencial(T[index_parcela], p[i])\
#                                   /potencialEquivalente(T[index_parcela], p[i], w[index_parcela])))
