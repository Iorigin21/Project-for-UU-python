from sympy import *
import math
#normalized equations

def electron_mobility_model(mu_n_max):
    return mu_n_max

def hole_mobility_model(mu_p_max):
    return mu_p_max

def Jn_m(phi_kp1,phi_k,n_kp1,n_k,mu_n,h_m, **kwargs):
    return mu_n/h_m*(phi_k-phi_kp1)*(n_kp1*exp(phi_k-phi_kp1)-n_k)/(exp(phi_k-phi_kp1)-1)

def Jp_m(phi_kp1,phi_k,p_kp1,p_k,mu_p,h_m, **kwargs):
    return mu_p/h_m*(phi_k-phi_kp1)*(p_k*exp(phi_k-phi_kp1)-p_kp1)/(exp(phi_k-phi_kp1)-1)

def Recombination(p_k,n_k,tau_n,tau_p):
    return (p_k*n_k-1)/(tau_p*(n_k+1)+tau_n*(p_k+1))

def F_poisson_equation(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k):
    return 2.0/(h_m+h_mm1)*((phi_kp1-phi_k)/h_m-(phi_k-phi_km1)/h_mm1)+p_k-n_k+dop_k

def F_electron_current_density(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k):
    return Jn_m(phi_kp1=phi_kp1,phi_k=phi_k,n_kp1=n_kp1,n_k=n_k,mu_n=mu_n,h_m=h_m)-Jn_m(phi_kp1=phi_k,phi_k=phi_km1,n_kp1=n_k,n_k=n_km1,mu_n=mu_n,h_m=h_mm1) - Recombination(p_k=p_k,n_k=n_k,tau_n=tau_n,tau_p=tau_p)

def F_hole_current_density(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k):
    return Jp_m(phi_kp1=phi_kp1,phi_k=phi_k,p_kp1=p_kp1,p_k=p_k,mu_p=mu_p,h_m=h_m)-Jp_m(phi_kp1=phi_k,phi_k=phi_km1,p_kp1=p_k,p_k=p_km1,mu_p=mu_p,h_m=h_mm1) + Recombination(p_k=p_k,n_k=n_k,tau_n=tau_n,tau_p=tau_p)


def Jn_m_eps(phi_kp1,phi_k,n_kp1,n_k,h_m,mu_n, **kwargs):
    return mu_n/h_m*(n_kp1*exp(phi_k-phi_kp1)-n_k)

def Jp_m_eps(phi_kp1,phi_k,p_kp1,p_k,h_m,mu_p, **kwargs):
    return mu_p/h_m*(p_k*exp(phi_k-phi_kp1)-p_kp1)
#for small phi_k-phi_kp1
def F_electron_current_density_eps(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k):
    return Jn_m_eps(phi_kp1=phi_kp1,phi_k=phi_k,n_kp1=n_kp1,n_k=n_k,mu_n=mu_n,h_m=h_m)-Jn_m_eps(phi_kp1=phi_k,phi_k=phi_km1,n_kp1=n_k,n_k=n_km1,mu_n=mu_n,h_m=h_mm1) - Recombination(p_k=p_k,n_k=n_k,tau_n=tau_n,tau_p=tau_p)

def F_hole_current_density_eps(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k):
    return Jp_m_eps(phi_kp1=phi_kp1,phi_k=phi_k,p_kp1=p_kp1,p_k=p_k,mu_p=mu_p,h_m=h_m)-Jp_m_eps(phi_kp1=phi_k,phi_k=phi_km1,p_kp1=p_k,p_k=p_km1,mu_p=mu_p,h_m=h_mm1) + Recombination(p_k=p_k,n_k=n_k,tau_n=tau_n,tau_p=tau_p)

##boundary condition
def BC_Phi_Va(Va,p_L):
    return Va - ln(p_L)

def BC_Neumann(N,Va):
    n   =   (N+sqrt(N**2 + 4))/2.0
    p   =   1.0/n
    phi =   Va - ln(p)
    return [n,p,phi]

#poisson equation
phi_kp1  = Symbol('phi_kp1')
phi_k    = Symbol('phi_k')
phi_km1  = Symbol('phi_km1')

h_m      = Symbol('h_m')
h_mm1    = Symbol('h_mm1')

n_k      = Symbol('n_k')
p_k      = Symbol('p_k')

dop_k      = Symbol('Dop_k')

n_kp1    = Symbol('n_kp1')
n_km1    = Symbol('n_km1')
p_kp1    = Symbol('p_kp1')
p_km1    = Symbol('p_km1')

mu_n     = Symbol('mu_n')
mu_p     = Symbol('mu_p')

tau_n    = Symbol('tau_n')
tau_p    = Symbol('tau_p')

def gen_jacobian_partial_diff_equations(equation,vars):
    if equation == "poisson":
        return diff(F_poisson_equation(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k),vars)
    elif equation == "electron":
        return diff(F_electron_current_density(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k),vars)
    elif equation == "hole":
        return diff(F_hole_current_density(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k),vars)
    elif equation == "hole_eps":
        return diff(F_hole_current_density_eps(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k),vars)
    elif equation == "electron_eps":
        return diff(F_electron_current_density_eps(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k),vars)
    elif equation == "Jn":
        return diff(Jn_m(phi_kp1=phi_kp1,phi_k=phi_k,n_kp1=n_kp1,n_k=n_k,h_m=h_m,mu_n=mu_n),vars)
    elif equation == "Jp":
        return diff(Jp_m(phi_kp1=phi_kp1,phi_k=phi_k,p_kp1=p_kp1,p_k=p_k,h_m=h_m,mu_p=mu_p),vars)
    elif equation == "Jn_eps":
        return diff(Jn_m_eps(phi_kp1=phi_kp1,phi_k=phi_k,n_kp1=n_kp1,n_k=n_k,h_m=h_m,mu_n=mu_n),vars)
    elif equation == "Jp_eps":
        return diff(Jp_m_eps(phi_kp1=phi_kp1,phi_k=phi_k,p_kp1=p_kp1,p_k=p_k,h_m=h_m,mu_p=mu_p), vars)
    elif equation == "Recombination":
        return diff(Recombination(p_k,n_k,tau_n,tau_p), vars)

#possion
poisson_phi_km1 = gen_jacobian_partial_diff_equations('poisson', phi_km1)
poisson_phi_k   = gen_jacobian_partial_diff_equations('poisson',phi_k)
poisson_phi_kp1 = gen_jacobian_partial_diff_equations('poisson',phi_kp1)
poisson_n_k = gen_jacobian_partial_diff_equations('poisson',n_k)
poisson_p_k = gen_jacobian_partial_diff_equations('poisson',p_k)
#electon
electron_phi_km1 = gen_jacobian_partial_diff_equations('electron', phi_km1)
electron_phi_k   = gen_jacobian_partial_diff_equations('electron',phi_k)
electron_phi_kp1 = gen_jacobian_partial_diff_equations('electron',phi_kp1)
electron_n_km1 = gen_jacobian_partial_diff_equations('electron', n_km1)
electron_n_k   = gen_jacobian_partial_diff_equations('electron', n_k)
electron_n_kp1 = gen_jacobian_partial_diff_equations('electron', n_kp1)
electron_p_k = gen_jacobian_partial_diff_equations('electron', p_k)
#hole
hole_phi_km1 = gen_jacobian_partial_diff_equations('hole', phi_km1)
hole_phi_k   = gen_jacobian_partial_diff_equations('hole', phi_k)
hole_phi_kp1 = gen_jacobian_partial_diff_equations('hole', phi_kp1)
hole_p_km1 = gen_jacobian_partial_diff_equations('hole', p_km1)
hole_p_k   = gen_jacobian_partial_diff_equations('hole', p_k)
hole_p_kp1 = gen_jacobian_partial_diff_equations('hole', p_kp1)
hole_n_k = gen_jacobian_partial_diff_equations('hole', n_k)


#electon
electron_phi_km1_eps = gen_jacobian_partial_diff_equations('electron_eps', phi_km1)
electron_phi_k_eps   = gen_jacobian_partial_diff_equations('electron_eps',phi_k)
electron_phi_kp1_eps = gen_jacobian_partial_diff_equations('electron_eps',phi_kp1)
electron_n_km1_eps = gen_jacobian_partial_diff_equations('electron_eps', n_km1)
electron_n_k_eps   = gen_jacobian_partial_diff_equations('electron_eps', n_k)
electron_n_kp1_eps = gen_jacobian_partial_diff_equations('electron_eps', n_kp1)
electron_p_k_eps = gen_jacobian_partial_diff_equations('electron_eps', p_k)
#hole
hole_phi_km1_eps = gen_jacobian_partial_diff_equations('hole_eps', phi_km1)
hole_phi_k_eps   = gen_jacobian_partial_diff_equations('hole_eps', phi_k)
hole_phi_kp1_eps = gen_jacobian_partial_diff_equations('hole_eps', phi_kp1)
hole_p_km1_eps = gen_jacobian_partial_diff_equations('hole_eps', p_km1)
hole_p_k_eps   = gen_jacobian_partial_diff_equations('hole_eps', p_k)
hole_p_kp1_eps = gen_jacobian_partial_diff_equations('hole_eps', p_kp1)
hole_n_k_eps = gen_jacobian_partial_diff_equations('hole_eps', n_k)


#Jn
Jn_phi_k   = gen_jacobian_partial_diff_equations('Jn', phi_k)
Jn_phi_kp1 = gen_jacobian_partial_diff_equations('Jn', phi_kp1)
Jn_n_k   = gen_jacobian_partial_diff_equations('Jn', n_k)
Jn_n_kp1 = gen_jacobian_partial_diff_equations('Jn', n_kp1)
#Jp
Jp_phi_k   = gen_jacobian_partial_diff_equations('Jp', phi_k)
Jp_phi_kp1 = gen_jacobian_partial_diff_equations('Jp', phi_kp1)
Jp_p_k   = gen_jacobian_partial_diff_equations('Jp', p_k)
Jp_p_kp1 = gen_jacobian_partial_diff_equations('Jp', p_kp1)


#Jn_eps
Jn_phi_k_eps   = gen_jacobian_partial_diff_equations('Jn_eps', phi_k)
Jn_phi_kp1_eps = gen_jacobian_partial_diff_equations('Jn_eps', phi_kp1)
Jn_n_k_eps   = gen_jacobian_partial_diff_equations('Jn_eps', n_k)
Jn_n_kp1_eps = gen_jacobian_partial_diff_equations('Jn_eps', n_kp1)
#Jp_eps
Jp_phi_k_eps   = gen_jacobian_partial_diff_equations('Jp_eps', phi_k)
Jp_phi_kp1_eps = gen_jacobian_partial_diff_equations('Jp_eps', phi_kp1)
Jp_p_k_eps   = gen_jacobian_partial_diff_equations('Jp_eps', p_k)
Jp_p_kp1_eps = gen_jacobian_partial_diff_equations('Jp_eps', p_kp1)


R_p_k  = gen_jacobian_partial_diff_equations('Recombination', p_k)
R_n_k  = gen_jacobian_partial_diff_equations('Recombination', n_k)

Print_Flag = True

def printFuncCode(equation,FunName):
    print('def %s(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k):' %(FunName))
    print('\treturn %s' % equation)
    print('')
    if Print_Flag:
        print('--------------%s公式如下----------------' % FunName)
        pprint(equation)
        print('----------------------------------------')

#gen for current
def printFuncCodeJn(equation,FunName):
    print('def %s(phi_kp1, phi_k, n_kp1, n_k,mu_n,h_m, **kwargs):' %(FunName))
    print('\treturn %s' % equation)
    print('')
    if Print_Flag:
        print('--------------%s公式如下----------------' % FunName)
        pprint(equation)
        print('----------------------------------------')

def printFuncCodeJp(equation,FunName):
    print('def %s(phi_kp1, phi_k, p_kp1, p_k,mu_p,h_m, **kwargs):' %(FunName))
    print('\treturn %s' % equation)
    print('')
    if Print_Flag:
        print('--------------%s公式如下----------------' % FunName)
        pprint(equation)
        print('----------------------------------------')

def printFuncCodeRecombination(equation, FunName):
    print('def %s(n_k, p_k,tau_n,tau_p, **kwargs):' % (FunName))
    print('\treturn %s' % equation)
    print('')
    if Print_Flag:
        print('--------------%s公式如下----------------' % FunName)
        pprint(equation)
        print('----------------------------------------')



def print_equations(equation):
    if not Print_Flag:
        return 0
    print('------------%s---------------' % equation)
    if equation == "poisson":
        pprint(F_poisson_equation(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k))
        print(latex(F_poisson_equation(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k)))
    elif equation == "electron":
        pprint(F_electron_current_density(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k))
    elif equation == "hole":
        pprint(F_hole_current_density(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k))
    elif equation == "hole_eps":
        pprint(F_hole_current_density_eps(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k))
    elif equation == "electron_eps":
        pprint(F_electron_current_density_eps(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k))
    elif equation == "Jn":
        pprint(Jn_m(phi_kp1=phi_kp1,phi_k=phi_k,n_kp1=n_kp1,n_k=n_k,h_m=h_m,mu_n=mu_n))
    elif equation == "Jp":
        pprint(Jp_m(phi_kp1=phi_kp1,phi_k=phi_k,p_kp1=p_kp1,p_k=p_k,h_m=h_m,mu_p=mu_p))
    elif equation == "Jn_eps":
        pprint(Jn_m_eps(phi_kp1=phi_kp1,phi_k=phi_k,n_kp1=n_kp1,n_k=n_k,h_m=h_m,mu_n=mu_n))
    elif equation == "Jp_eps":
        pprint(Jp_m_eps(phi_kp1=phi_kp1,phi_k=phi_k,p_kp1=p_kp1,p_k=p_k,h_m=h_m,mu_p=mu_p))
    elif equation == "Recombination":
        pprint(Recombination(p_k=p_k,n_k=n_k,tau_n=tau_n,tau_p=tau_p))

if __name__ == "__main__":
    printFuncCode(poisson_phi_km1,'poisson_phi_km1')
    printFuncCode(poisson_phi_k,'poisson_phi_k')
    printFuncCode(poisson_phi_kp1,'poisson_phi_kp1')
    printFuncCode(poisson_n_k,'poisson_n_k')
    printFuncCode(poisson_p_k,'poisson_p_k')

    printFuncCode(electron_phi_km1,'electron_phi_km1')
    printFuncCode(electron_phi_k,'electron_phi_k')
    printFuncCode(electron_phi_kp1,'electron_phi_kp1')
    printFuncCode(electron_n_km1,'electron_n_km1')
    printFuncCode(electron_n_k,'electron_n_k')
    printFuncCode(electron_n_kp1,'electron_n_kp1')
    printFuncCode(electron_p_k,'electron_p_k')

    printFuncCode(hole_phi_km1,'hole_phi_km1')
    printFuncCode(hole_phi_k,'hole_phi_k')
    printFuncCode(hole_phi_kp1,'hole_phi_kp1')
    printFuncCode(hole_p_km1,'hole_p_km1')
    printFuncCode(hole_p_k,'hole_p_k')
    printFuncCode(hole_p_kp1,'hole_p_kp1')
    printFuncCode(hole_n_k,'hole_n_k')

    printFuncCode(electron_phi_km1_eps,'electron_phi_km1_eps')
    printFuncCode(electron_phi_k_eps,'electron_phi_k_eps')
    printFuncCode(electron_phi_kp1_eps,'electron_phi_kp1_eps')
    printFuncCode(electron_n_km1_eps,'electron_n_km1_eps')
    printFuncCode(electron_n_k_eps,'electron_n_k_eps')
    printFuncCode(electron_n_kp1_eps,'electron_n_kp1_eps')
    printFuncCode(electron_p_k_eps,'electron_p_k_eps')

    printFuncCode(hole_phi_km1_eps,'hole_phi_km1_eps')
    printFuncCode(hole_phi_k_eps,'hole_phi_k_eps')
    printFuncCode(hole_phi_kp1_eps,'hole_phi_kp1_eps')
    printFuncCode(hole_p_km1_eps,'hole_p_km1_eps')
    printFuncCode(hole_p_k_eps,'hole_p_k_eps')
    printFuncCode(hole_p_kp1_eps,'hole_p_kp1_eps')
    printFuncCode(hole_n_k_eps,'hole_n_k_eps')

    #Jn
    printFuncCodeJn(Jn_phi_k, 'Jn_phi_k')
    printFuncCodeJn(Jn_phi_kp1, 'Jn_phi_kp1')
    printFuncCodeJn(Jn_n_k, 'Jn_n_k')
    printFuncCodeJn(Jn_n_kp1, 'Jn_n_kp1')

    #Jp
    printFuncCodeJp(Jp_phi_k, 'Jp_phi_k')
    printFuncCodeJp(Jp_phi_kp1, 'Jp_phi_kp1')
    printFuncCodeJp(Jp_p_k, 'Jp_p_k')
    printFuncCodeJp(Jp_p_kp1, 'Jp_p_kp1')

    printFuncCodeJn(Jn_phi_k_eps, 'Jn_phi_k_eps')
    printFuncCodeJn(Jn_phi_kp1_eps, 'Jn_phi_kp1_eps')
    printFuncCodeJn(Jn_n_k_eps, 'Jn_n_k_eps')
    printFuncCodeJn(Jn_n_kp1_eps, 'Jn_n_kp1_eps')

    printFuncCodeJp(Jp_phi_k_eps, 'Jp_phi_k_eps')
    printFuncCodeJp(Jp_phi_kp1_eps, 'Jp_phi_kp1_eps')
    printFuncCodeJp(Jp_p_k_eps, 'Jp_p_k_eps')
    printFuncCodeJp(Jp_p_kp1_eps, 'Jp_p_kp1_eps')

    printFuncCodeRecombination(R_p_k, 'R_p_k')
    printFuncCodeRecombination(R_n_k, 'R_n_k')


    print_equations('poisson')
    print_equations('electron')
    print_equations('electron_eps')
    print_equations('hole')
    print_equations('hole_eps')
    print_equations('Jn')
    print_equations('Jn_eps')
    print_equations('Jp')
    print_equations('Jp_eps')
    print_equations('Recombination')