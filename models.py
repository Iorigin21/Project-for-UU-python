from math import *


def BC_Neumann(N,Va):
    n   =   (N+sqrt(N**2 + 4))/2.0
    p   =   1.0/n
    phi =   Va - log(p)
    return [n,p,phi]

def electron_mobility_model(mu_n_max):
    return mu_n_max

def hole_mobility_model(mu_p_max):
    return mu_p_max

def Jn_m(phi_kp1,phi_k,n_kp1,n_k,mu_n,h_m, **kwargs):
    return mu_n/h_m*(phi_k-phi_kp1)*(n_kp1*exp(phi_k-phi_kp1)-n_k)/(exp(phi_k-phi_kp1)-1)

def Jp_m(phi_kp1,phi_k,p_kp1,p_k,mu_p,h_m, **kwargs):
    return mu_p/h_m*(phi_k-phi_kp1)*(p_k*exp(phi_k-phi_kp1)-p_kp1)/(exp(phi_k-phi_kp1)-1)

def Recombination(p_k,n_k,tau_n,tau_p, **kwargs):
    return (p_k*n_k-1)/(tau_n*(n_k+1)+tau_p*(p_k+1))

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




'''######################Auto generated Code By models_sympy.py##########################'''
def poisson_phi_km1(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k):
	return 2.0/(h_mm1*(h_m + h_mm1))

def poisson_phi_k(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k):
	return 2.0*(-1/h_mm1 - 1/h_m)/(h_m + h_mm1)

def poisson_phi_kp1(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k):
	return 2.0/(h_m*(h_m + h_mm1))

def poisson_n_k(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k):
	return -1

def poisson_p_k(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k):
	return 1

def electron_phi_km1(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k):
	return -mu_n*n_k*(-phi_k + phi_km1)*exp(-phi_k + phi_km1)/(h_mm1*(exp(-phi_k + phi_km1) - 1)) + mu_n*(-phi_k + phi_km1)*(n_k*exp(-phi_k + phi_km1) - n_km1)*exp(-phi_k + phi_km1)/(h_mm1*(exp(-phi_k + phi_km1) - 1)**2) - mu_n*(n_k*exp(-phi_k + phi_km1) - n_km1)/(h_mm1*(exp(-phi_k + phi_km1) - 1))

def electron_phi_k(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k):
	return mu_n*n_k*(-phi_k + phi_km1)*exp(-phi_k + phi_km1)/(h_mm1*(exp(-phi_k + phi_km1) - 1)) - mu_n*(-phi_k + phi_km1)*(n_k*exp(-phi_k + phi_km1) - n_km1)*exp(-phi_k + phi_km1)/(h_mm1*(exp(-phi_k + phi_km1) - 1)**2) + mu_n*(n_k*exp(-phi_k + phi_km1) - n_km1)/(h_mm1*(exp(-phi_k + phi_km1) - 1)) + mu_n*n_kp1*(phi_k - phi_kp1)*exp(phi_k - phi_kp1)/(h_m*(exp(phi_k - phi_kp1) - 1)) - mu_n*(-n_k + n_kp1*exp(phi_k - phi_kp1))*(phi_k - phi_kp1)*exp(phi_k - phi_kp1)/(h_m*(exp(phi_k - phi_kp1) - 1)**2) + mu_n*(-n_k + n_kp1*exp(phi_k - phi_kp1))/(h_m*(exp(phi_k - phi_kp1) - 1))

def electron_phi_kp1(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k):
	return -mu_n*n_kp1*(phi_k - phi_kp1)*exp(phi_k - phi_kp1)/(h_m*(exp(phi_k - phi_kp1) - 1)) + mu_n*(-n_k + n_kp1*exp(phi_k - phi_kp1))*(phi_k - phi_kp1)*exp(phi_k - phi_kp1)/(h_m*(exp(phi_k - phi_kp1) - 1)**2) - mu_n*(-n_k + n_kp1*exp(phi_k - phi_kp1))/(h_m*(exp(phi_k - phi_kp1) - 1))

def electron_n_km1(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k):
	return mu_n*(-phi_k + phi_km1)/(h_mm1*(exp(-phi_k + phi_km1) - 1))

def electron_n_k(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k):
	return -p_k/(tau_n*(p_k + 1) + tau_p*(n_k + 1)) + tau_p*(n_k*p_k - 1)/(tau_n*(p_k + 1) + tau_p*(n_k + 1))**2 - mu_n*(-phi_k + phi_km1)*exp(-phi_k + phi_km1)/(h_mm1*(exp(-phi_k + phi_km1) - 1)) - mu_n*(phi_k - phi_kp1)/(h_m*(exp(phi_k - phi_kp1) - 1))

def electron_n_kp1(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k):
	return mu_n*(phi_k - phi_kp1)*exp(phi_k - phi_kp1)/(h_m*(exp(phi_k - phi_kp1) - 1))

def electron_p_k(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k):
	return -n_k/(tau_n*(p_k + 1) + tau_p*(n_k + 1)) + tau_n*(n_k*p_k - 1)/(tau_n*(p_k + 1) + tau_p*(n_k + 1))**2

def hole_phi_km1(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k):
	return -mu_p*p_km1*(-phi_k + phi_km1)*exp(-phi_k + phi_km1)/(h_mm1*(exp(-phi_k + phi_km1) - 1)) + mu_p*(-p_k + p_km1*exp(-phi_k + phi_km1))*(-phi_k + phi_km1)*exp(-phi_k + phi_km1)/(h_mm1*(exp(-phi_k + phi_km1) - 1)**2) - mu_p*(-p_k + p_km1*exp(-phi_k + phi_km1))/(h_mm1*(exp(-phi_k + phi_km1) - 1))

def hole_phi_k(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k):
	return mu_p*p_km1*(-phi_k + phi_km1)*exp(-phi_k + phi_km1)/(h_mm1*(exp(-phi_k + phi_km1) - 1)) - mu_p*(-p_k + p_km1*exp(-phi_k + phi_km1))*(-phi_k + phi_km1)*exp(-phi_k + phi_km1)/(h_mm1*(exp(-phi_k + phi_km1) - 1)**2) + mu_p*(-p_k + p_km1*exp(-phi_k + phi_km1))/(h_mm1*(exp(-phi_k + phi_km1) - 1)) + mu_p*p_k*(phi_k - phi_kp1)*exp(phi_k - phi_kp1)/(h_m*(exp(phi_k - phi_kp1) - 1)) - mu_p*(phi_k - phi_kp1)*(p_k*exp(phi_k - phi_kp1) - p_kp1)*exp(phi_k - phi_kp1)/(h_m*(exp(phi_k - phi_kp1) - 1)**2) + mu_p*(p_k*exp(phi_k - phi_kp1) - p_kp1)/(h_m*(exp(phi_k - phi_kp1) - 1))

def hole_phi_kp1(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k):
	return -mu_p*p_k*(phi_k - phi_kp1)*exp(phi_k - phi_kp1)/(h_m*(exp(phi_k - phi_kp1) - 1)) + mu_p*(phi_k - phi_kp1)*(p_k*exp(phi_k - phi_kp1) - p_kp1)*exp(phi_k - phi_kp1)/(h_m*(exp(phi_k - phi_kp1) - 1)**2) - mu_p*(p_k*exp(phi_k - phi_kp1) - p_kp1)/(h_m*(exp(phi_k - phi_kp1) - 1))

def hole_p_km1(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k):
	return -mu_p*(-phi_k + phi_km1)*exp(-phi_k + phi_km1)/(h_mm1*(exp(-phi_k + phi_km1) - 1))

def hole_p_k(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k):
	return n_k/(tau_n*(p_k + 1) + tau_p*(n_k + 1)) - tau_n*(n_k*p_k - 1)/(tau_n*(p_k + 1) + tau_p*(n_k + 1))**2 + mu_p*(-phi_k + phi_km1)/(h_mm1*(exp(-phi_k + phi_km1) - 1)) + mu_p*(phi_k - phi_kp1)*exp(phi_k - phi_kp1)/(h_m*(exp(phi_k - phi_kp1) - 1))

def hole_p_kp1(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k):
	return -mu_p*(phi_k - phi_kp1)/(h_m*(exp(phi_k - phi_kp1) - 1))

def hole_n_k(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k):
	return p_k/(tau_n*(p_k + 1) + tau_p*(n_k + 1)) - tau_p*(n_k*p_k - 1)/(tau_n*(p_k + 1) + tau_p*(n_k + 1))**2

def electron_phi_km1_eps(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k):
	return -mu_n*n_k*exp(-phi_k + phi_km1)/h_mm1

def electron_phi_k_eps(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k):
	return mu_n*n_k*exp(-phi_k + phi_km1)/h_mm1 + mu_n*n_kp1*exp(phi_k - phi_kp1)/h_m

def electron_phi_kp1_eps(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k):
	return -mu_n*n_kp1*exp(phi_k - phi_kp1)/h_m

def electron_n_km1_eps(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k):
	return mu_n/h_mm1

def electron_n_k_eps(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k):
	return -p_k/(tau_n*(p_k + 1) + tau_p*(n_k + 1)) + tau_p*(n_k*p_k - 1)/(tau_n*(p_k + 1) + tau_p*(n_k + 1))**2 - mu_n*exp(-phi_k + phi_km1)/h_mm1 - mu_n/h_m

def electron_n_kp1_eps(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k):
	return mu_n*exp(phi_k - phi_kp1)/h_m

def electron_p_k_eps(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k):
	return -n_k/(tau_n*(p_k + 1) + tau_p*(n_k + 1)) + tau_n*(n_k*p_k - 1)/(tau_n*(p_k + 1) + tau_p*(n_k + 1))**2

def hole_phi_km1_eps(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k):
	return -mu_p*p_km1*exp(-phi_k + phi_km1)/h_mm1

def hole_phi_k_eps(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k):
	return mu_p*p_km1*exp(-phi_k + phi_km1)/h_mm1 + mu_p*p_k*exp(phi_k - phi_kp1)/h_m

def hole_phi_kp1_eps(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k):
	return -mu_p*p_k*exp(phi_k - phi_kp1)/h_m

def hole_p_km1_eps(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k):
	return -mu_p*exp(-phi_k + phi_km1)/h_mm1

def hole_p_k_eps(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k):
	return n_k/(tau_n*(p_k + 1) + tau_p*(n_k + 1)) - tau_n*(n_k*p_k - 1)/(tau_n*(p_k + 1) + tau_p*(n_k + 1))**2 + mu_p/h_mm1 + mu_p*exp(phi_k - phi_kp1)/h_m

def hole_p_kp1_eps(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k):
	return -mu_p/h_m

def hole_n_k_eps(phi_kp1, phi_k, phi_km1, n_kp1, n_k, n_km1, p_kp1,p_k,p_km1,h_m, h_mm1,mu_n,mu_p,tau_n,tau_p,dop_k):
	return p_k/(tau_n*(p_k + 1) + tau_p*(n_k + 1)) - tau_p*(n_k*p_k - 1)/(tau_n*(p_k + 1) + tau_p*(n_k + 1))**2

def Jn_phi_k(phi_kp1, phi_k, n_kp1, n_k,mu_n,h_m, **kwargs):
	return mu_n*n_kp1*(phi_k - phi_kp1)*exp(phi_k - phi_kp1)/(h_m*(exp(phi_k - phi_kp1) - 1)) - mu_n*(-n_k + n_kp1*exp(phi_k - phi_kp1))*(phi_k - phi_kp1)*exp(phi_k - phi_kp1)/(h_m*(exp(phi_k - phi_kp1) - 1)**2) + mu_n*(-n_k + n_kp1*exp(phi_k - phi_kp1))/(h_m*(exp(phi_k - phi_kp1) - 1))

def Jn_phi_kp1(phi_kp1, phi_k, n_kp1, n_k,mu_n,h_m, **kwargs):
	return -mu_n*n_kp1*(phi_k - phi_kp1)*exp(phi_k - phi_kp1)/(h_m*(exp(phi_k - phi_kp1) - 1)) + mu_n*(-n_k + n_kp1*exp(phi_k - phi_kp1))*(phi_k - phi_kp1)*exp(phi_k - phi_kp1)/(h_m*(exp(phi_k - phi_kp1) - 1)**2) - mu_n*(-n_k + n_kp1*exp(phi_k - phi_kp1))/(h_m*(exp(phi_k - phi_kp1) - 1))

def Jn_n_k(phi_kp1, phi_k, n_kp1, n_k,mu_n,h_m, **kwargs):
	return -mu_n*(phi_k - phi_kp1)/(h_m*(exp(phi_k - phi_kp1) - 1))

def Jn_n_kp1(phi_kp1, phi_k, n_kp1, n_k,mu_n,h_m, **kwargs):
	return mu_n*(phi_k - phi_kp1)*exp(phi_k - phi_kp1)/(h_m*(exp(phi_k - phi_kp1) - 1))

def Jp_phi_k(phi_kp1, phi_k, p_kp1, p_k,mu_p,h_m, **kwargs):
	return mu_p*p_k*(phi_k - phi_kp1)*exp(phi_k - phi_kp1)/(h_m*(exp(phi_k - phi_kp1) - 1)) - mu_p*(phi_k - phi_kp1)*(p_k*exp(phi_k - phi_kp1) - p_kp1)*exp(phi_k - phi_kp1)/(h_m*(exp(phi_k - phi_kp1) - 1)**2) + mu_p*(p_k*exp(phi_k - phi_kp1) - p_kp1)/(h_m*(exp(phi_k - phi_kp1) - 1))

def Jp_phi_kp1(phi_kp1, phi_k, p_kp1, p_k,mu_p,h_m, **kwargs):
	return -mu_p*p_k*(phi_k - phi_kp1)*exp(phi_k - phi_kp1)/(h_m*(exp(phi_k - phi_kp1) - 1)) + mu_p*(phi_k - phi_kp1)*(p_k*exp(phi_k - phi_kp1) - p_kp1)*exp(phi_k - phi_kp1)/(h_m*(exp(phi_k - phi_kp1) - 1)**2) - mu_p*(p_k*exp(phi_k - phi_kp1) - p_kp1)/(h_m*(exp(phi_k - phi_kp1) - 1))

def Jp_p_k(phi_kp1, phi_k, p_kp1, p_k,mu_p,h_m, **kwargs):
	return mu_p*(phi_k - phi_kp1)*exp(phi_k - phi_kp1)/(h_m*(exp(phi_k - phi_kp1) - 1))

def Jp_p_kp1(phi_kp1, phi_k, p_kp1, p_k,mu_p,h_m, **kwargs):
	return -mu_p*(phi_k - phi_kp1)/(h_m*(exp(phi_k - phi_kp1) - 1))

def Jn_phi_k_eps(phi_kp1, phi_k, n_kp1, n_k,mu_n,h_m, **kwargs):
	return mu_n*n_kp1*exp(phi_k - phi_kp1)/h_m

def Jn_phi_kp1_eps(phi_kp1, phi_k, n_kp1, n_k,mu_n,h_m, **kwargs):
	return -mu_n*n_kp1*exp(phi_k - phi_kp1)/h_m

def Jn_n_k_eps(phi_kp1, phi_k, n_kp1, n_k,mu_n,h_m, **kwargs):
	return -mu_n/h_m

def Jn_n_kp1_eps(phi_kp1, phi_k, n_kp1, n_k,mu_n,h_m, **kwargs):
	return mu_n*exp(phi_k - phi_kp1)/h_m

def Jp_phi_k_eps(phi_kp1, phi_k, p_kp1, p_k,mu_p,h_m, **kwargs):
	return mu_p*p_k*exp(phi_k - phi_kp1)/h_m

def Jp_phi_kp1_eps(phi_kp1, phi_k, p_kp1, p_k,mu_p,h_m, **kwargs):
	return -mu_p*p_k*exp(phi_k - phi_kp1)/h_m

def Jp_p_k_eps(phi_kp1, phi_k, p_kp1, p_k,mu_p,h_m, **kwargs):
	return mu_p*exp(phi_k - phi_kp1)/h_m

def Jp_p_kp1_eps(phi_kp1, phi_k, p_kp1, p_k,mu_p,h_m, **kwargs):
	return -mu_p/h_m

def R_p_k(n_k, p_k,tau_n,tau_p, **kwargs):
	return n_k/(tau_n*(p_k + 1) + tau_p*(n_k + 1)) - tau_n*(n_k*p_k - 1)/(tau_n*(p_k + 1) + tau_p*(n_k + 1))**2

def R_n_k(n_k, p_k,tau_n,tau_p, **kwargs):
	return p_k/(tau_n*(p_k + 1) + tau_p*(n_k + 1)) - tau_p*(n_k*p_k - 1)/(tau_n*(p_k + 1) + tau_p*(n_k + 1))**2

