from cell1d import Cell1D
import matplotlib.pyplot as plt
import numpy as np

class SEMI_DEVICE():
    def __init__(self,mesh_data,*args,**kwargs):
        self.mesh_point     =   mesh_data['mesh_points']
        self.doping_data    =   mesh_data['doping_data']
        self.cell_array     =   []
        self.num_cell       =   None

        self.plot           = plt.figure()
        #initial structure
        self.gen_meshed_structure()

    def gen_meshed_structure(self):
        len_mesh = len(self.mesh_point)
        assert len_mesh > 1
        for i in range(len_mesh-1):
            x0 = self.mesh_point[i]
            x1 = self.mesh_point[i+1]
            doping0 =   self.doping_data[i]
            doping1 =   self.doping_data[i+1]
            if i == 0:
                cell0 = Cell1D(X_kp1=x1,X_k=x0,Dop_kp1 = doping1,Dop_k=doping0,border= "min")
            elif i == len_mesh - 2:
                cell0 = Cell1D(X_kp1=x1,X_k=x0,Dop_kp1 = doping1,Dop_k=doping0,border= "max")
            else:
                cell0 = Cell1D(X_kp1=x1, X_k=x0, Dop_kp1=doping1, Dop_k=doping0, border=None)
            self.cell_array.append(cell0)
        self.num_cell = len(self.cell_array)

    def draw_structure(self):
        plt.figure(1)
        plt.plot(self.mesh_point,self.doping_data,'ro-')
        plt.show()

    def traversal_cells(self):
        for cell0 in self.cell_array:
            print('Cell%3d: x0(%.2e) Dop_k0:%.2e phi0:%.2e(%.2e) n0:%.2e(%.2e) p0:%.2e(%.2e)' %(cell0.id,cell0.X_k,cell0.Dop_k,cell0.Phi_k,cell0.phi_k,cell0.N_k,cell0.n_k,cell0.P_k,cell0.p_k) )
            print('-------: x1(%.2e) Dop_k1:%.2e phi1:%.2e(%.2e) n1:%.2e(%.2e) p1:%.2e(%.2e)' % (
            cell0.X_kp1, cell0.Dop_k, cell0.Phi_kp1, cell0.phi_kp1, cell0.N_kp1, cell0.n_kp1, cell0.P_kp1,
            cell0.p_kp1))

    def set_material_parameters(self,parameters):
        for cell0 in self.cell_array:
            cell0.initialize_parameters(parameters)

    def set_electrode_voltage(self,vol):
        self.cell_array[-1].set_electrode_voltage(vol[1])
        self.cell_array[0].set_electrode_voltage(vol[0])

    def update(self,T0,normalize_operation = False):
        for (index,cell0) in enumerate(self.cell_array):
            #normalization
            cell0.set_variables(T0[2,index],T0[0,index],T0[1,index],T0[2,index+1],T0[0,index+1],T0[1,index+1],normalize_operation)
            cell0.comput_diffs()

    def get_initial_variable(self):
        num_of_cell = len(self.cell_array)
        T0 = np.zeros((3,num_of_cell+1))
        cell_last = None
        for (index,cell0) in enumerate(self.cell_array):
            T0[2,index] = cell0.phi0
            T0[0,index] = cell0.n0
            T0[1,index] = cell0.p0
            cell_last = cell0
        T0[2, -1] = cell_last.phi1
        T0[0, -1] = cell_last.n1
        T0[1, -1] = cell_last.p1
        return T0

    def current(self):
        return self.cell_array[0].current_density

def plot_device(device,plt,i,label):
    X = []
    Y = []
    Y2 = []
    Y3 = []
    Y4 = []
    Y5 = []
    Y6 = []
    cell_array = device.cell_array
    for cell0 in cell_array:
        X.append(cell0.X_k)
        Y.append(cell0.Phi_k)
        Y2.append(cell0.E)
        Y3.append(cell0.N_k)
        Y4.append(cell0.P_k)
        Y5.append(cell0.electron_current_density)
        Y6.append(cell0.hole_current_density)
    XE = X.copy()
    X.append(cell_array[-1].X_kp1)
    X = np.array(X)
    XE = np.array(XE)
    X = X*1e4   #transform to unit um
    XE = XE*1e4 #transform to unit um
    Y.append(cell_array[-1].Phi_kp1)
    Y3.append(cell_array[-1].N_kp1)
    Y4.append(cell_array[-1].P_kp1)
    color_list = ['red','orange','olive','green','darkslategray','aqua','slategrey','darkblue','indigo','purple']
    ax1 = plt.subplot(3, 2, 1)
    ax1.plot(X, Y, linestyle = '-',color = color_list[i % len(color_list)],label = label )
    ax2 = plt.subplot(3, 2, 2)
    ax2.plot(XE, Y2, linestyle = '-',color = color_list[i % len(color_list)],label = label )
    ax3 = plt.subplot(3, 2, 3)
    plt.yscale('log')
    ax3.plot(X, Y3, linestyle = '-',color = color_list[i % len(color_list)],label = label )
    ax4 = plt.subplot(3, 2, 4)
    plt.yscale('log')
    ax4.plot(X, Y4, linestyle = '-',color = color_list[i % len(color_list)],label = label )
    ax5 = plt.subplot(3, 2, 5)
    ax5.plot(XE, Y5, linestyle='-', color=color_list[i % len(color_list)], label=label)
    ax6 = plt.subplot(3, 2, 6)
    ax6.plot(XE, Y6, linestyle='-', color=color_list[i % len(color_list)], label=label)

    ax1.set_xlabel(r'Location, X (${\mu}m$)')
    ax2.set_xlabel(r'Location, X (${\mu}m$)')
    ax3.set_xlabel(r'Location, X (${\mu}m$)')
    ax4.set_xlabel(r'Location, X (${\mu}m$)')
    ax5.set_xlabel(r'Location, X (${\mu}m$)')
    ax6.set_xlabel(r'Location, X (${\mu}m$)')

    ax1.set_ylabel(r'Potential, $\varphi$ (V)')
    ax2.set_ylabel(r'ElectriField, $E$ (V/cm)')
    ax3.set_ylabel(r'Electron Density, $n$ ($cm^{-3}$)')
    ax4.set_ylabel(r'Hole Density, $p$ ($cm^{-3}$)')
    ax5.set_ylabel(r'Electron Current Density, $J_{n}$ ($A/cm^{2}$)')
    ax6.set_ylabel(r'Hole Current Density, $J_{p}$ ($A/cm^{2}$)')
    ax1.set_title('ElectroStaticPotential')
    ax2.set_title('ElectricField')
    ax3.set_title('Electron Density')
    ax4.set_title('Hole Density')
    ax5.set_title('Electron Current Density')
    ax6.set_title('Hole Current Density')
    #plt.subplots_adjust(left=None, bottom=None, right=None, top=0.5, wspace=0.5, hspace=0.5)

    ax1.legend()
    ax2.legend()
    ax3.legend()
    ax4.legend()
    ax5.legend()
    ax6.legend()

def plot_energyband(device,plt,title):
    x = []
    Fermi_e = []
    Fermi_h = []
    Ec_list = list()
    Ev_list = list()
    cell_array = device.cell_array
    for cell0 in cell_array:
        x.append(cell0.X_k)
        Fermi_e.append(cell0.Fermi_electron)
        Fermi_h.append(cell0.Fermi_hole)
        Ec_list.append(cell0.Ec)
        Ev_list.append(cell0.Ev)
    xe = x.copy()
    xe = np.array(xe)*1e4 #transform to unit:cm
    color_list = ['red', 'orange', 'olive', 'green', 'darkslategray', 'aqua', 'slategrey', 'darkblue', 'indigo',
                  'purple']
    ax1 = plt.gca()
    ax1.plot(xe, Fermi_e, linestyle='--', color=color_list[1], label=r"$E_{Fn}$")
    ax1.plot(xe, Fermi_h, linestyle='--', color=color_list[2], label=r"$E_{Fp}$")
    ax1.plot(xe, Ec_list, linestyle='-', color=color_list[3], label=r"$E_{C}$")
    ax1.plot(xe, Ev_list, linestyle='-', color=color_list[4], label=r"$E_{V}$")
    # ax1.suptitle(title)
    ax1.set_xlabel(r'Location, $X$ (${\mu}m$)')
    ax1.set_ylabel(r'EnergyBand, $E$ (eV)')
    ax1.legend()

def print_welcome_text():
    print('''
*******************************************************************************
***                          1D PN diode Device                             ***
***                          Version R-2022.03v1.0                          ***
***                        (x86_64, Python 3.10)                            ***
***                                                                         ***
***                       Copyright (C) For Everyone                        ***
***                                                                         ***
***     The project for Advanced Scientific Programming in python           ***
***     Help to understand the working process of TCAD                      ***
***                                                                         ***
*******************************************************************************
    ''')

if __name__ == "__main__":
    L1  = 1e-4 #cm
    L0  = 0
    L_junction = 0.5e-4
    Dop_P      = 1e17
    Dop_N      = 1e18
    mesh_points = np.linspace(L0,L1,30)
    mesh_data  = np.empty(len(mesh_points))
    for (index,point) in enumerate(mesh_points):
        if point> L_junction:
            mesh_data[index] = Dop_N
        else:
            mesh_data[index] = -1*Dop_P
    mesh_data  = {'mesh_points':mesh_points,'doping_data':mesh_data}
    pn_diode = SEMI_DEVICE(mesh_data)
    pn_diode.draw_structure()
    pn_diode.traversal_cells()