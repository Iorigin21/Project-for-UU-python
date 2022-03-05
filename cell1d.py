from normalize import *
from models import *

EPS =  1e-10

class Cell1D():
    cell_count = 0
    def __init__(self, X_kp1, X_k, Dop_kp1,Dop_k,border = None ):
        self.x0 = X_k / LOC_NORMAL
        self.x1 = X_kp1 / LOC_NORMAL
        # already normalized
        self.hm0 = (self.x1 - self.x0)

        self.dop0   = Dop_k/DOPPING_NORMAL
        self.dop1   = Dop_kp1/DOPPING_NORMAL

        #count the num of the cells
        Cell1D.cell_count += 1
        self.id = Cell1D.cell_count

        #border
        self.border = border

        #initialize the variables
        self.initialize_variables()
        self.initialize_BC()

    def initialize_variables(self, initial_values = None):
        #初始条件设置
        if self.dop0 > 0:
            self.n0 = self.dop0
            self.p0 = 1/self.dop0
            # 初始条件
            self.phi0 = log(self.n0)

        else:
            self.n0 = -1/self.dop0
            self.p0 = -1*self.dop0
            self.phi0 = - log(self.p0)

        if self.dop1 > 0:
            self.n1 = self.dop1
            self.p1 = 1/self.dop1
            self.phi1 = log(self.n1)
        else:
            self.n1 = -1/self.dop1
            self.p1 = -1*self.dop1
            self.phi1 = - log(self.p1)


    def initialize_BC(self):
        # Va applied voltage
        self.v0 = 0/PHI_NORMAL

    def initialize_parameters(self,paramters):
        #normalized
        self.par_mu_n_max = paramters['mu_n_max']/MOBILITY_NROMAL
        self.par_mu_p_max = paramters['mu_p_max']/MOBILITY_NROMAL

        self.par_tau_n_max = paramters['tau_n_max'] / TIME_NORMAL
        self.par_tau_p_max = paramters['tau_p_max'] / TIME_NORMAL

        self.mu_n = electron_mobility_model(self.par_mu_n_max)
        self.mu_p = hole_mobility_model(self.par_mu_p_max)

        self.tau_n = electron_mobility_model(self.par_tau_n_max)
        self.tau_p = hole_mobility_model(self.par_tau_p_max)

    def comput_diffs(self):
        par_dict = self.normalized_values()
        if abs(self.phi0-self.phi1)<EPS:
            # print("EPS")
            self.Jn = Jn_m_eps(**par_dict)
            self.Jn_phi_k = Jn_phi_k_eps(**par_dict)
            self.Jn_phi_kp1 = Jn_phi_kp1_eps(**par_dict)
            self.Jn_n_k = Jn_n_k_eps(**par_dict)
            self.Jn_n_kp1 = Jn_n_kp1_eps(**par_dict)


            self.Jp = Jp_m_eps(**par_dict)
            self.Jp_phi_k = Jp_phi_k_eps(**par_dict)
            self.Jp_phi_kp1 = Jp_phi_kp1_eps(**par_dict)
            self.Jp_p_k = Jp_p_k_eps(**par_dict)
            self.Jp_p_kp1 = Jp_p_kp1_eps(**par_dict)

        else:
            self.Jn = Jn_m(**par_dict)
            self.Jn_phi_k = Jn_phi_k(**par_dict)
            self.Jn_phi_kp1 = Jn_phi_kp1(**par_dict)
            self.Jn_n_k = Jn_n_k(**par_dict)
            self.Jn_n_kp1 = Jn_n_kp1(**par_dict)

            self.Jp = Jp_m(**par_dict)
            self.Jp_phi_k = Jp_phi_k(**par_dict)
            self.Jp_phi_kp1 = Jp_phi_kp1(**par_dict)
            self.Jp_p_k = Jp_p_k(**par_dict)
            self.Jp_p_kp1 = Jp_p_kp1(**par_dict)

        self.R_n_k = R_n_k(**par_dict)
        self.R_p_k = R_p_k(**par_dict)

        self.R     = Recombination(**par_dict)

        #print(self.Jn_phi_k)
    @property
    def X_k(self):
        return self.x0 * LOC_NORMAL

    @property
    def X_kp1(self):
        return self.x1 * LOC_NORMAL

    @property
    def hm(self):
        return self.hm0

    @property
    def Phi_k(self):
        return self.phi0*PHI_NORMAL

    @property
    def phi_kp1(self):
        return self.phi1

    @property
    def phi_k(self):
        return self.phi0

    @property
    def Phi_kp1(self):
        return self.phi1 * PHI_NORMAL

    @property
    def N_k(self):
        return self.n0*CARRIER_NORMAL

    @property
    def n_k(self):
        return self.n0

    @property
    def N_kp1(self):
        return self.n1 * CARRIER_NORMAL

    @property
    def n_kp1(self):
        return self.n1

    @property
    def P_k(self):
        return self.p0*CARRIER_NORMAL

    @property
    def p_k(self):
        return self.p0

    @property
    def Fermi_electron(self):
        return (log(self.n0)/log(e) + self.phi0) * PHI_NORMAL

    @property
    def Fermi_hole(self):
        return -1.0*(log(self.p0) / log(e) - self.phi0) * PHI_NORMAL


    @property
    def Ec(self):
        return self.phi0 * PHI_NORMAL - phi_ref + x_aff + Eg

    @property
    def Ev(self):
        return self.phi0 * PHI_NORMAL - phi_ref + x_aff

    @property
    def P_kp1(self):
        return self.p1 * CARRIER_NORMAL

    @property
    def p_kp1(self):
        return self.p1

    @property
    def Dop_k(self):
        return self.dop0 * DOPPING_NORMAL

    @property
    def dop_k(self):
        return self.dop0

    @property
    def dop_kp1(self):
        return self.dop1

    @property
    def Dop_kp1(self):
        return self.dop1 * DOPPING_NORMAL

    @property
    def electron_current_density(self):
        self.comput_diffs()
        return -1*self.Jn*CURRENT_DEN_NORMAL

    @property
    def hole_current_density(self):
        self.comput_diffs()
        return -1*self.Jp*CURRENT_DEN_NORMAL

    @property
    def current_density(self):
        self.comput_diffs()
        return -1*(self.Jn+self.Jp) * CURRENT_DEN_NORMAL

    @property
    def E(self):
        return (self.phi0-self.phi1)/self.hm*E_NORMAL

    @property
    def va(self):
        return self.v0

    @property
    def Va(self):
        return self.v0 * PHI_NORMAL

    def set_electrode_voltage(self,Va):
        #设置边界条件
        self.v0 = Va/PHI_NORMAL
        if self.border == "min":
            (self.n0, self.p0, self.phi0) = BC_Neumann(self.dop_k, self.v0)
        elif self.border == "max":
            (self.n1, self.p1, self.phi1) = BC_Neumann(self.dop_kp1, self.v0)
            # print('SET:Max Boundary Condition V= %.3f V' % Va)


    def set_variables(self,phi0,n0,p0,phi1,n1,p1,normalize_operation=False):
        if normalize_operation:
            #归一化以后的
            #如果本结构是左边界，则不赋值
            if self.border != "min":
                self.phi0 = phi0/PHI_NORMAL
                self.n0   = n0/CARRIER_NORMAL
                self.p0   = p0/CARRIER_NORMAL
            # 如果本结构是右边界，则不赋值
            if self.border != "max":
                self.phi1 = phi1/PHI_NORMAL
                self.n1   = n1/CARRIER_NORMAL
                self.p1   = p1/CARRIER_NORMAL
        else:
            #如果本结构是左边界，则不赋值
            if self.border != "min":
                self.phi0 = phi0
                self.n0   = n0
                self.p0   = p0
            # 如果本结构是右边界，则不赋值
            if self.border != "max":
                self.phi1 = phi1
                self.n1   = n1
                self.p1   = p1

    def normalized_values(self):
        normalized_values_dict = {
            'phi_kp1'   :   self.phi1,
            'phi_k'     :   self.phi0,

            'n_kp1'     :   self.n1,
            'n_k'       :   self.n0,

            'p_kp1'     :   self.p1,
            'p_k'       :   self.p0,

            'x_kp1'     :   self.x1,
            'x_k'       :   self.x0,
            'h_m'       :   self.hm0,

            'dop_k'     :   self.dop0,
            'dop_kp1'   :   self.dop1,

            'mu_n'      :   self.mu_n,
            'mu_p'      :   self.mu_p,

            'tau_n': self.tau_n,
            'tau_p': self.tau_p
        }
        return normalized_values_dict


