from sympy import *
import sys
import matplotlib.pyplot as plt
class Memoize:
    def __init__(self, f):
        self.f = f
        self.memo = {}
    def __call__(self, *args):
        if not args in self.memo:
            self.memo[args] = self.f(*args)
        #Warning: You may wish to do a deepcopy here if returning objects
        return self.memo[args]
sys.setrecursionlimit(1000)


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
    partial_three = sum([sum([D(i,j,V,p0)*C(n-i,m+1-j,V,p0) for j in range(1,i+2)]) for i in range(1,n+1)])
    return -1/(2*D(0,1,V,p0)) *(partial_one+ partial_two + 2*partial_three)
@Memoize
def D(n,m,V,p0):
    if m >= n+2:
        return 0
    if m<=0:
        return 0
    if n == 0 and m == 1:
        return -sqrt(2*W(0,2,V,p0))
    
    partial_one = -2*W(2*n,2*m,V,p0)
    partial_two = (2*m+1)*D(n,m+1,V,p0)
    partial_three = sum([sum([D(i,j,V,p0)*D(n-i,m+1-j,V,p0) for j in range(1,i+2)]) for i in range(1,n)])
    partial_four = sum([sum([C(i,j,V,p0)*C(n-i-1,m-j,V,p0) for j in range(0,i+2)]) for i in range(0,n)])
    return -1/(2*D(0,1,V,p0))*(partial_one+partial_two+ partial_three+partial_four) 
def Enminusone(n,V,p0,var,value):

    partial_one = -D(n,1,V,p0)+2*W(2*n,0,V,p0)
    partial_two = sum([C(i,0,V,p0)*C(n-i-1,0,V,p0) for i in range(0,n)])
    return (n-1,(1/2 * (partial_one- partial_two).subs([(var,value)])*3**(-n)))
def get_energy_order(n,option):

    #ls =  3**(4/3)/(2**(8/3)) * 3
    #ls = -2/9
    V,p0,Eminustwo,var,value = potential(option)
    ls = 3*Eminustwo.subs([(var,value)])
    for i in range(0,n):
        ls+= Enminusone(i,V,p0,var,value)[1]
    return ls

def plot(n,option):
    V,p0,Eminustwo,var,value = potential(option)
    ls = 3*Eminustwo.subs([(var,value)])
    n_list = [i for i in range(n)]
    Elist = [0 for i in range(n)]
    for j in range(n):
        
        if j == 0:
            Elist[j] = ls
        else:
            Elist[j] = Elist[j-1] + Enminusone(j-1,V,p0,var,value)[1]
        print(Elist[j])
        print(j)
    return(n_list,Elist)
def potential(option):
    p = Symbol('p',positive =True)
    var = Symbol('var',positive =True)
    if option == 'Coulomb':
        p0 = 1/(4*var)
        potential = -var/p
        value = 1/3**(3/2)
    elif option == 'Quartic':
        potential = var* p**4
        p0 = 1/(2**(2/3)*var**(1/6))
        value = 3
    elif option == 'Linear':
        potential = var*p
        p0 =1/(var**(1/3)*2**(2/3))
        value =  2**(7/2)/3**(1/2)
    elif option == 'Quintic':
        potential = 3**(3/2)*p**5
        p0 = 1/(2**(2/7)*3**(3/14)*5**(1/7))
        value = 1
    elif option == 'Harmonic':
        potential =var *p**2
        p0 = 1/(2**(3/4)*var**(1/4))
        value = 1/2

    V = (1/(8*p**2) + potential)
    Eminustwo = V.subs(p,p0)
    return (V,p0,Eminustwo,var,value)

#print(get_energy_order(20,'Coulomb'))
print(plot(50,'Linear'))