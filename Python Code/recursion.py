from sympy import *
import sys
import matplotlib.pyplot as plt
import functools
import collections
import cProfile
import re
import numpy as np


class Memoize(dict):
    def __init__(self, func):
         self.func = func 
    def __call__(self, *args):
        return self[args]
 
    def __missing__(self, key):
        result = self[key] = self.func(*key)
        return result


@Memoize
def W(n,m,V,p0):
    p = Symbol('p',positive =True)

    if n == 0:
        if m == 0:
            return -1/(2*p0**2)
        if m == 2:
            return diff(V,p,2).subs(p,p0)/factorial(2)
        else:
            return 0
    if n== 1:
        if m == 1:
            return 1/(p0**3)
        if m == 3:
            return diff(V,p,3).subs(p,p0)/factorial(3)
        else:
            return 0
        
    if m == n-2:
        return (-1)**n *3*(n-1)/(8*p0**n)
    elif m == n:
        return (-1)**(n+1)*(n+1)/(2*p0**(n+2))
    elif m == n+2:
        return diff(V,p,n+2).subs(p,p0)/factorial(n+2)
    else:
        return 0
@Memoize
def C(n,m,V,p0):
    
    if m >= n+2:
        return 0
    if m<0:
        return 0
    partial_one = -2 *W(2*n+1,2*m+1,V,p0)
    partial_two = 2*(m+1)*C(n,m+1,V,p0)
    partial_three = np.sum([np.sum([D(i,j,V,p0)*C(n-i,m+1-j,V,p0) for j in range(1,i+2)]) for i in range(1,n+1)])
    result = -1/(2*D(0,1,V,p0)) *(partial_one+ partial_two + 2*partial_three)
    return result
@Memoize
def D(n,m,V,p0):
    if m >= n+2:
        return 0
    if m<=0:
        return 0
    if n == 0 and m == 1:
        a = str(-sqrt(2*W(0,2,V,p0)))
        return -sqrt(2*W(0,2,V,p0))
    
    partial_one = -2*W(2*n,2*m,V,p0)
    partial_two = (2*m+1)*D(n,m+1,V,p0)
    partial_three = np.sum([np.sum([D(i,j,V,p0)*D(n-i,m+1-j,V,p0) for j in range(1,i+2)]) for i in range(1,n)])
    partial_four = np.sum([np.sum([C(i,j,V,p0)*C(n-i-1,m-j,V,p0) for j in range(0,i+2)]) for i in range(0,n)])
    result = -1/(2*D(0,1,V,p0))*(partial_one+partial_two+ partial_three+partial_four) 
    return result
def Enminusone(n,V,p0,var,value,k):

    partial_one = -D(n,1,V,p0)+2*W(2*n,0,V,p0)
    partial_two = np.sum([C(i,0,V,p0)*C(n-i-1,0,V,p0) for i in range(0,n)])
    return (n-1,(1/2 * N((partial_one- partial_two).subs([(var,value)])*k**(-n))))
def get_energy_order(n,option):

    V,p0,Eminustwo,var,value = potential(option)
    ls = 3*Eminustwo.subs([(var,value)])
    for i in range(0,n):
        ls+= Enminusone(i,V,p0,var,value)[1]
    return ls

def plot(n,option):
    V,p0,Eminustwo,var,value,k = potential(option)
    ls = Eminustwo.subs([(var,value)])
    n_list = [i for i in range(n)]
    Elist = [0 for i in range(n)]
    for j in range(n):
        
        if j == 0:
            Elist[j] = ls
        else:
            Elist[j] = collect(Elist[j-1] + Enminusone(j-1,V,p0,var,value,k)[1],k)
        print(j)
    return(n_list,[2*N(Elist[i]) for i in range(len(Elist))])
def potential(option):
    k = Symbol('k',positive = True)
    p = Symbol('p',positive = True )
    var = Symbol('var',positive = True)
    if option == 'Coulomb':
        p0 =1/(4*var)
        potential = - var/p
        value = 1/3**(3/2)
    elif option == 'Quartic':
        potential = var* p**4
        p0 = 1/(2**(2/3)*var**(1/6))
        value = N
    elif option == 'Linear':
        potential = var*p
        p0 =1/((var)**(1/3)*2**(2/3))
        value =  2**(7/2)/3**(1/2)
    elif option == 'Quintic':
        potential =var*p**5
        p0 = 1/(2**(2/7)*(var)**(1/7)*5**(1/7))
        value =  3**(3/2)
    elif option == 'Harmonic':
        potential =-var *p**2
        p0 = 1/(2**(3/4)*var**(1/4))
        value = 1/2
    elif option == 'Type 1':
        potential = -var * p**(-1/5)
        p0 = (5/(4*var))**(5/9)
        value = 2**(1.7)/(5*sqrt(5)**(1/5))
    elif option == 'Type 2':
        potential = -var * p**(-4/5)
        p0 = (5/(16*var))**(5/6)
        value = 2**(0.8)/(3*sqrt(3)**(4/5))
    elif option == 'Sqrt':
        potential = var*p**(1/2)
        p0 = (1/(2*var))**(2/5)
        value = 1/(2*3**(3/4))
    elif option == 'Type 3':
        potential = - var*p**(-1.5)
        p0 = 1/(36*var**2)
        value = 1/(2*3**(7/4))
    elif option == 'Log':
        potential = var*log(3**(1/2)*p)
        p0 = 1/(4*var)**(1/2)
        value = 1/(3*2)
    elif option == 'Inverse Sqrt':
        potential = -var*p**(-1/2)
        p0 = 1/((2*var)**(2/3))
        value = 1/3**(5/4)
    elif option == 'Charm':
        potential = -3**(-3/2)/(2*p) +p*3**(-1/2)/(2)
        p0  =0.8372665140289375792830123
        value =1
    elif option == 'Hulten':
        a = 0.1
        potential = -0.1/7*exp(-7**(1/2)*0.1*p)/(1-exp(-0.1*7**(1/2)*p))
        p0 =5.5106
        value = 1
    elif option == 'Morse':
        potential = 10/3*(exp(-2*3**(1/2)*p)-2*exp(-3**(1/2)*p))
        p0 = 0.443191
        value =1
    elif option == 'Martin':
        potential = 6.8698/(3*2)*(p*3**(1/2))**(0.1)-8.064/(3*2)
        p0 = 1.41299
        value =1
    elif option == 'Dressed Coulomb':
        potential = -1/9 * 1/(9*p**2+1)**(1/2)
        p0 =6.77452771823243
        value = 1
    elif option == 'Spiked Harmonic':
        potential = p**2/2 + 0.01/(p*3**(3/2)*2)
        p0 = 0.707588
        value = 1
    elif option == 'Max Min':
        potential = p**2/2 -3**(1/2)*p**3/(10*2)
        p0 = 0.746252
        value = 1

    k = 3
    V = (1/(8*p**2) + potential)
    Eminustwo = V.subs(p,p0)*k
    return (V,p0,Eminustwo,var,value,k)
#V,p0,Eminustwo,vars,value = potential('Coulomb')
#print(Enminusone(10,V,p0,vars,value))
a = plot(10,"Max Min")
print(a)
#cProfile.run('re.compile(print(plot(50,"Coulomb")))')
