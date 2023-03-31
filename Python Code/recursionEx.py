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
def a(n,V,p0):
    partial_one = -sum([a(k,V,p0)*T(n-k,0,V,p0) for k in range(1,n)])
    return 1/(T(0,0,V,p0))*(2*C(n-1,0,V,p0)+ partial_one)
@Memoize
def S(n,m,V,p0):
    partial_one = sum([a(k,V,p0)*T(n-k+1,m+1,V,p0) for k in range(1,n-m+1)])
    return partial_one -2*C(n,m+1,V,p0)
@Memoize
def T(n,m,V,p0):
    partial_one = sum([a(k,V,p0)*S(n-k,m,V,p0) for k in range(1,n-m+1)])
    return partial_one - 2*D(n,m+1,V,p0)
@Memoize
def C(n,m,V,p0):
    
    if m >= n+2:
        return 0
    if m<0:
        return 0
    partial_one = -S(n,m,V,p0)-2 *W(2*n+1,2*m+1,V,p0)
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
        return -sqrt(2*W(0,2,V,p0))
    
    partial_one = -T(n,m,V,p0)-2*W(2*n,2*m,V,p0)
    partial_two = (2*m+1)*D(n,m+1,V,p0)
    partial_three = np.sum([np.sum([D(i,j,V,p0)*D(n-i,m+1-j,V,p0) for j in range(1,i+2)]) for i in range(1,n)])
    partial_four = np.sum([np.sum([C(i,j,V,p0)*C(n-i-1,m-j,V,p0) for j in range(0,i+2)]) for i in range(0,n)])
    result = -1/(2*D(0,1,V,p0))*(partial_one+partial_two+ partial_three+partial_four) 
    return result
def Enminusone(n,V,p0,var,value,k):

    partial_one =T(n,0,V,p0) -D(n,1,V,p0)+2*W(2*n,0,V,p0)
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
    return(n_list,[N(Elist[i]) for i in range(len(Elist))])
def potential(option):
    k = Symbol('k',positive = True)
    p = Symbol('p',positive = True )
    var = Symbol('var',positive = True)
    if option == 'Coulomb':
        p0 =1/(4*var)
        potential = - var/p
        value = 1/3**(3/2)
    if option == 'Sqrt':
        potential = var*p**(1/2)
        p0 = (1/(2*var))**(2/5)
        value = 1/(3**(3/4))
    elif option == 'Type 3':
        potential = - var*p**(-1.5)
        p0 = 1/(36*var**2)
        value = 1/(3**(7/4))
    elif option == 'Type 4':
        potential = var*p**(3/20)
        p0 = (5/(3*var))**(20/43)
        value = 1/3**(37/40)
    elif option == 'Log':
        potential = var*log(3**(1/2)*p)
        p0 = 1/(4*var)**(1/2)
        value = 1/(3)
    elif option == 'Quartic':
        potential = var* p**4
        p0 = 1/(2**(2/3)*var**(1/6))
        value = 3
    elif option == 'Quintic':
        potential =var*p**5
        p0 = 1/(2**(2/7)*(var)**(1/7)*5**(1/7))
        value =  3**(3/2)
    elif option == 'Linear':
        potential = var*p
        p0 =1/((var)**(1/3)*2**(2/3))
        value =  2**(7/2)/3**(1/2)
    elif option == 'Dressed Coulomb':
        k = 3
        lambd = 1
        potential = -1/k * 1/(k*p**2+lambd**2)**(1/2)
        p0 = nsolve(diff(1/(8*p**2) + potential,p),p,1)
        value = 1
    elif option == 'Cornell':
        val = 0.92
        a = 0.52
        b = 1/(2.34**2)
        k = 5
        potential =val* (-a*k**(-3/2)/(p) +b*p*k**(-1/2))
        p0 = nsolve(diff(1/(8*p**2) + potential,p),p,1)
        value =1
    elif option == 'Spiked Harmonic':
        a = 0.000001
        k=3
        potential = p**2 + a/(p**(4)*k**(3))
        p0 = nsolve(diff(1/(8*p**2) + potential,p),p,1)
        value = 1
    k = 3
    V = (1/(8*p**2) + potential)
    Eminustwo = V.subs(p,p0)*k
    return (V,p0,Eminustwo,var,value,k)
#V,p0,Eminustwo,vars,value = potential('Coulomb')
#print(Enminusone(10,V,p0,vars,value))
print(plot(10,"Spiked Harmonic"))
#cProfile.run('re.compile(print(plot(50,"Coulomb")))')
