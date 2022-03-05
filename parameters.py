import math

# material settings
MATERIAL = "Silicon"


#common const
epsilon0    =   8.854178e-14    # [F/cm]
T           =   300.0           # [K]
q           =   1.6e-19         # [C]
##Boltzmann constant
k           =   1.380649e-23    # [J/K]
Vt          =   k*T/q           # [V]


########Material Parameter
#if MATERIAL == "Silicon":
#intrisic density
x_aff       =   4.05 #[eV]
Eg0         =   1.12 #[eV]
Eg          =   Eg0 #[eV] #no bandgap narrowing effect
Nc          =   2.89e19 #[cm^-3]
Nv          =   3.14e19 #[cm^-3]
ni          =   (Nc*Nv)**(0.5)*math.exp(-1.0*Eg/2.0/Vt) #[cm^-3]
epsilon     =   12.0    #[1]
#mobility
mu_n_max    =   1330    #[cm^2/V/s]
mu_p_max    =   495    #[cm^2/V/s]
## diffusion coefficient(Einstein relation)
D_n         =   mu_n_max*Vt
D_p         =   mu_p_max*Vt
D_0         =   1*Vt
##life time of the carriers
tau_n_max   =   0.5e-6
tau_p_max   =   0.5e-6 #0.5e-10

##potential reference
phi_ref     = x_aff + Eg/2.0

material_parameters = {
    'mu_n_max':mu_n_max,
    'mu_p_max':mu_p_max,
    'tau_n_max':tau_n_max,
    'tau_p_max':tau_p_max
}