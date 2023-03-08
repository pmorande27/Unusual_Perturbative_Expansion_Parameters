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
def W(n,m):
    p = Symbol('p',positive =True)
    Abar = Symbol('Abar',positive=True)
    ebar = Symbol('ebar',positive = True)
    #potential = 3* p**4
   
    potential = -ebar/p
    

    V = 1/(8*p**2) + potential
    p0 =1/(4*ebar)
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
def C(n,m):
    if m >= n+2:
        return 0
    if m<0:
        return 0
    partial_one = -2 *W(2*n+1,2*m+1)
    partial_two = 2*(m+1)*C(n,m+1)
    partial_three =0
    for i in range(1,n+1):
        for j in range(1,i+2):
            partial_three += D(i,j)*C(n-i,m+1-j)
    return simplify(-1/(2*D(0,1)) *(partial_one+ partial_two + 2*partial_three))
@Memoize
def D(n,m):
    #print('a')
    if m >= n+2:
        return 0
    if m<=0:
        return 0
    if n == 0 and m == 1:
        return -sqrt(2*W(0,2))
    
    partial_one = -2*W(2*n,2*m)
    partial_two = (2*m+1)*D(n,m+1)
    partial_three = 0
    for i in range(1,n):
        for j in range(1,i+2):
            partial_three += D(i,j)*D(n-i,m+1-j)
    partial_four = 0
    for i in range(0,n):
        for j in range(0,i+2):
            a = n-i-1
            b = m-j
            partial_four += C(i,j)*C(a,b)
    return -1/(2*D(0,1))*(partial_one+partial_two+ partial_three+partial_four) 
def Enminusone(n):
    ebar = Symbol('ebar',positive = True)
    e = Symbol('e',positive = True)
    k = Symbol('k',positive = True)

    partial_one = -D(n,1)+2*W(2*n,0)
    partial_two = 0
    for i in range(0,n):
    
        partial_two += C(i,0)*C(n-i-1,0)
    return (n-1,(1/2 * (partial_one- partial_two).subs(ebar,1/(3)**(3/2))*3**(-n)))
def get_energy_order(n):
    ls =  -2/9
    for i in range(0,n):
        ls+= Enminusone(i)[1]
    return ls

def plot(n):
    n_list = [i for i in range(n)]
    Elist = [0 for i in range(n)]
    for j in range(n):
        print(j)
        Elist[j] = get_energy_order(j)
    return(n_list,Elist)


print(plot(20))