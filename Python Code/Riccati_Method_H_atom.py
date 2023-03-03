from sympy import *
import matplotlib.pyplot as plt
# Setting Up all variables and functions of use
x = Symbol('x',positive = True)
ebar = Symbol('ebar',positive = True)
r = Symbol('r',positive = True)
Abar = Symbol('Abar')
e = Symbol('e',positive = True)
V = Abar* r**4
#V = -ebar**2/r
k = Symbol('k', positive = True)
A = Symbol('A',positive=True)
Veff = 1/(8* r**2) + V
r0 = 4**(1/(-5))* (Abar *4)**(1/(-6))
#r0 = 1/(4*ebar**2)
#print(Veff)
Veffbar = factor(simplify(Veff.subs(r,x+r0)- Veff.subs(r,r0)))
#print(Veffbar)


#Getting the first coefficients
Eminustwo = Veff.subs(r,r0)
#print(Eminustwo)
phiminusone =  simplify(-sqrt(2* Veffbar))
#print(phiminusone)
Eminusone = -1/2* diff(phiminusone,x).subs(x,0) - 1/(2* r0**2)
#print(Eminusone)
phizero = simplify((1/2 * diff(phiminusone,x) + 1/(2*(x+r0)**2) + Eminusone)/(-phiminusone))
Ezero = -1/2 *diff(phizero,x).subs(x,0) - 1/2* phizero.subs(x,0)**2 + 3/(8* r0**2)
phione = simplify((1/2* diff(phizero,x) + 1/2* phizero**2 - 3/(8* (x+r0)**2) + Ezero)/(-phiminusone))
Eone = -1/2* diff(phione,x).subs(x,0) - 1/2 *(phione.subs(x,0)* phizero.subs(x,0) + phione.subs(x,0)* phizero.subs(x,0))
#print(diff(phione,x).subs(x,0))
phis = [phiminusone,phizero,phione]
phis0 = [phizero,phione]
Elist = [Eminustwo,Eminusone,Ezero,Eone]
#print(Elist)
#print(phis)


# Dynamic Programming approach
def approximate_to_order(Elistinit,phiszero,order):
    iter = order -4
    """if order < 4:
        return Elistinit[0:order]"""
    Elist =Elistinit.copy()
    phis0 = phiszero.copy()
    n = 2
    for i in range(iter):
        Rel_E = Elist[len(Elist)-1]
        sum_term = 0
        for m in range(0,n):
            sum_term += 1/2 * phis0[m]*phis0[n-m-1]
        prime_term = simplify(1/2* diff(phis0[len(phis0)-1],x))
        phi_next = simplify((Rel_E + sum_term + prime_term)/(-phiminusone))
        #print(phi_next)
        phis0.append(phi_next.copy())

        E_prime_term = -1/2 *diff(phis0[len(phis0)-1],x).subs(x,0)
        sum_term_E = 0
        for j in range(0,n+1):
            sum_term_E+= -1/2* phis0[j].subs(x,0)*phis0[n-j].subs(x,0)
        E_next = sum_term_E+ E_prime_term
        Elist.append(E_next.copy())
        n+= 1
    klist = [k**(-i) for i in range(-2,n)]
    if order < 4:
        length = len(Elist) - (4-order)
    else:
        length = len(Elist)
    Energy_list = [klist[i]*Elist[i] for i in range(length)]
    Energy_value = sum([Energy_list[i].subs([(Abar,1/k),(k,3)]) for i in range(len(Energy_list))])
    #Energy_value = sum([Energy_list[i].subs([(ebar,e/k),(k,3)]) for i in range(len(Energy_list))])
    #Real_value = -e**4/2
    Real_value =2.393644
    error = abs((Energy_value-Real_value)/Real_value)
    accuracy = 1-error
    #print(accuracy)
    return (Energy_list,error)

#print(approximate_to_order(Elist,phiszero=phis0,order= 5))
def plot(n):
    n_list = [i for i in range(1,n+1)]
    error_list = [0]*n
    for j in range(n):
        #print(j+1)
        error_list[j] = approximate_to_order(Elistinit=Elist,phiszero=phis0,order=n_list[j])[1]
    print(n_list,error_list)

plot(10)
