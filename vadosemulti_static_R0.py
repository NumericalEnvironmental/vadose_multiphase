############################################################################################################################################
#
# VadoseMulti.py: a python script for simulating 2D multiphase (air/water/LNAPL) redistribution in a soil column
#
# * employ method of lines, with S_w and S_NAPL, as the primary variables in the coupled PDEs
# * pressure is calculated as a separate dependent variable by solving separate pressure equations for unsaturated and saturated conditions
# * effective Van Genuchten alpha parameter is scaled to the proportions of liquids in each cell
# * set NAPL residual saturation to zero, except for systems that are initially NAPL-soaked
# * capillary pressure not considered for liquid-saturated cells (only buoyancy forces and relative permeability)
#
# - methods are generally set up for broadcasting
#
############################################################################################################################################

from numpy import *
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from copy import *

# constants ...
g = 9.807           # gravitational acceleration
epsilon = 1e-6      # offset to prevent divide-by-zero error in some routines
min_Sn = 1e-10      # minimum NAPL saturation for calculation of relative permeability
S_crit_0 = 0.999    # critical liquid saturation above which pressure will be calculated via linear extrapolation
S_crit_f =1.00001   # critical liquid supersaturation above which pressure will be calculated by storage compression (>= 1.0 + dS_0[i] for all i)
min_Ss = 0.0001     # maximum specific storage value

###########################################################
#
# classes
#
###########################################################


class Liquid:
    def __init__(self,name,rho,u,beta):
        self.name = name
        self.rho = rho              # density
        self.u = u                  # viscosity
        self.beta = beta            # capillary scaling factors (2-element vector for NAPL-water and NAPL-air)

class Cell:
    def __init__(self,z,poros,k,alpha,n,rs,dS_0,dS_f,Ss_0,Ss_f,vol):
        self.z = z                  # elevation
        self.poros = poros          # porosity
        self.k = k                  # intrinsic permeability
        self.alpha = alpha          # Van Genuchten parameters
        self.n = n
        self.m = 1.0 - 1.0/n
        self.rs = rs                # residual saturation; two-element vector (water and NAPL)      
        self.dS_0 = dS_0            # minimum (0) and maximum (f) fluid oversaturations (>1.0) for defining specific storage curve
        self.dS_f = dS_f
        self.Ss_0 = Ss_0            # specific storage curve values corresponding to dS_0 and dS_f; a logarithmic progression is assumed
        self.Ss_f = Ss_f
        self.vol = vol              # cell volume; set to a large number to model a fixed pressure boundary
    def P_sat(self,S_w,S_n,S_t,liquids):
        # compute saturated fluid pressures (i.e., P > 0.) implied by interpolation of specific storage between two oversaturation endpoints
        dS = (S_t - 1.0)*(S_t>1.0) + self.dS_0*(S_t<=1.0)
        x = (log10(dS)-log10(self.dS_0))/(log10(self.dS_f)-log10(self.dS_0))
        Ss_p = 10**(log10(self.Ss_0) + x*(log10(self.Ss_f)-log10(self.Ss_0)))
        Ss = Ss_p*(Ss_p >= min_Ss) + min_Ss * (min_Ss<0.0001)
        return (dS * self.poros * (S_w*liquids.rho[0]/S_t + S_n*liquids.rho[1]/S_t) * g/Ss) * (S_t>1.0) + 0.0
    def P_unsat(self,S_t,alpha_eff):
        # return capillary pressure (as a negative number) for total liquid saturation
        St_0 = S_t*(S_t<=S_crit_0) + S_crit_0*(S_t>S_crit_0)
        return -exp(log(exp(-log(St_0)/self.m) - 1.0)/self.n-log(alpha_eff))*(S_t<=S_crit_0) + 0.0
    def P_bridge(self,P_crit_0,P_crit_f,S_t):
        # return pressure (linearly interpolated) spanning transistion from unsaturated to saturated
        return (P_crit_0+(P_crit_f-P_crit_0)*(S_t-S_crit_0)/(S_crit_f-S_crit_0))*((S_t > S_crit_0).astype(int)*(S_t < S_crit_f).astype(int)) + 0.0
    def Pressure(self,S_w,S_n,liquids):
        # return fluid pressure corresponding to saturation
        # note: S_t and S_w + S_n can be used interchangeably for bridged and saturated conditions with insignificant differences
        S_e,S_t = self.SatEff(S_w,S_n)
        alpha_eff = self.AlphaEff(S_w,S_n,1.0-S_w-S_n,liquids)
        P_crit_0 = self.P_unsat(S_crit_0,alpha_eff)         # pressure value at saturation value where linear extrapolation is employed
        P_crit_f = self.P_sat(S_w,S_n,S_crit_f,liquids)
        P_unsatd = self.P_unsat(S_t,alpha_eff)
        P_bridged = self.P_bridge(P_crit_0,P_crit_f,S_t)
        P_satd = self.P_sat(S_w,S_n,S_t,liquids)
        return P_unsatd + P_bridged + P_satd
    def SatEff(self,S_w,S_n):
        # effective liquid and total liquids saturation vectors
        S_e = array([(S_w - self.rs.T[0])/(1.0 - self.rs.T[0]),(S_n - self.rs.T[1])/(1.0 - self.rs.T[1])]).T
        S_t = (S_w + S_n - self.rs.T[0] - self.rs.T[1])/(1.0 - self.rs.T[0] - self.rs.T[1])
        return S_e,S_t
    def PermEff(self,S_w,S_n,liquids):
        # return effective permeabilities (relative permeability x intrinsic permeability) for both fluids (as a vector)
        Se_0,St_0 = self.SatEff(S_w,S_n)                      # update effective fluid saturations
        Se_w = Se_0.T[0] * (Se_0.T[0] <= 1.0) + 1.0 * (Se_0.T[0] > 1.0)
        Se_n = Se_0.T[1] * (Se_0.T[1] <= 1.0) + 1.0 * (Se_0.T[1] > 1.0)
        S_t = St_0 * (St_0 <= 1.0) + 1.0 * (St_0 > 1.0)
        k0 = Se_w**0.5 * (1.0 - (1.0 - Se_w**(1.0/self.m))**self.m)**2
        k1 = (S_t - Se_w)**0.5 * ((1.0-Se_w**(1.0/self.m))**self.m - (1.0-S_t**(1.0/self.m))**self.m)**2
        kr = array([k0/liquids.u[0],(k1*(Se_n>min_Sn) + 0.0)/liquids.u[1]])     # note: min_Sn limit used to suppress spurious NAPL movement from numerical residuals
        return kr * self.k
    def AlphaEff(self,sw0,sn0,sa0,liquids):
        # function to return weighted Van Genuchten alpha factor, based on phase saturations
        sw = 0.0*(sw0<0.0) + sw0*(sw0>=0.0)*(sw0<=1.0) + 1.0*(sw0>1.0)
        sn = 0.0*(sn0<0.0) + sn0*(sn0>=0.0)*(sn0<=1.0) + 1.0*(sn0>1.0)
        sa = 0.0*(sa0<0.0) + sa0*(sa0>=0.0)*(sa0<=1.0) + 1.0*(sa0>1.0)        
        aw = sw + sa
        nw = sw + sn
        na = sn + sa
        aw0 = 1.0/abs(1.0 - aw + epsilon)
        nw0 = 1.0/abs(1.0 - nw + epsilon)
        na0 = 1.0/abs(1.0 - na + epsilon)
        norm = aw0 + nw0 + na0
        alpha_eff = aw0/norm*self.alpha + nw0/norm*liquids.beta[1,0]*self.alpha + na0/norm*liquids.beta[1,1]*self.alpha
        return alpha_eff
    def Q(self,phase,S_w_0,S_n_0,A,dz,liquids):                 # this function returns a vector of net fluid fluxes of fluid 'phase' for all cell connections
        # suppress unphysical low saturations (<rs)
        S_w = S_w_0*(S_w_0>=self.rs.T[0]) + self.rs.T[0]*(S_w_0<self.rs.T[0])
        S_n = S_n_0*(S_n_0>=self.rs.T[1]) + self.rs.T[1]*(S_n_0<self.rs.T[1])
        # fluid fluxes (n_cells - 1 members; positive number = downward-directed flow)
        P = self.Pressure(S_w,S_n,liquids)
        grad = ((P[:-1]-P[1:]) + liquids.rho[phase]*g*(self.z[:-1]-self.z[1:]))/dz
        upstream_K = (0.5+0.5*sign(grad))*self.PermEff(S_w,S_n,liquids)[phase][:-1] + (0.5-0.5*sign(grad))*self.PermEff(S_w,S_n,liquids)[phase][1:]
        flux = A * upstream_K * grad
        # append 0's to ends of flux array (to faciliate processing by CoupledODES method)
        flux = append(array([0.0]),flux)
        flux = append(flux,array([0.0]))
        return flux
    def CoupledODEs(self,y,t,liquids,n_cells,A,dz):                              # systems of coupled ODEs
        dy = zeros(len(y),float)
        # water volumetric balance equations for all cells;
        dy[:n_cells] = 1.0/(self.vol*self.poros) \
            * (self.Q(0,y[:n_cells],y[n_cells:],A,dz,liquids)[:-1] - self.Q(0,y[:n_cells],y[n_cells:],A,dz,liquids)[1:])
        # NAPL volumetric balance equations for all cells;
        dy[n_cells:] = 1.0/(self.vol*self.poros) \
            * (self.Q(1,y[:n_cells],y[n_cells:],A,dz,liquids)[:-1] - self.Q(1,y[:n_cells],y[n_cells:],A,dz,liquids)[1:])
        return dy
    def ODEHandler(self,t,S_w,S_n,liquids,n_cells,A,dz):
        # set up system of ODEs (lump S_w and S_n into 'y') and call ODE solver routine
        y = zeros(2*n_cells,float)  # 1st portion = S_w, 2nd portion = S_n
        for i in xrange(n_cells):
            y[i] = S_w[i]
            y[i+n_cells] = S_n[i]
        y_t = odeint(self.CoupledODEs,y,t,mxstep=5000,args=(liquids,n_cells,A,dz))
        return y_t


###########################################################
#
# misc. support functions
#
###########################################################


def ProcessInput(prefix):
    # read input file
    name = []
    rho = zeros(2,float)
    u = zeros(2,float)
    beta = zeros((2,2),float)
    input_file = open(prefix + '_input.txt','r')
    flag_read = [0,0,0]
    i = 0
    for line in input_file:
        line_input = line.split()
        if line_input and line_input[0] == 'liquids':
            flag_read = [1,0,0]
            i = 0
        if line_input and line_input[0] == 'column':
            flag_read = [0,1,0]
            i = 0        
        if line_input and line_input[0] == 'cells':
            flag_read = [0,0,1]
            i = 0
        # set up liquid objects
        if flag_read[0] and i > 0 and i <= 2:
            name.append(line_input[0])
            rho[i-1] = float(line_input[1])
            u[i-1] = float(line_input[2])
            beta[i-1] = [float(line_input[3]),float(line_input[4])]
        # read column properties
        if flag_read[1] and i > 0 and i <= 3:
            if line_input[0] == 'span':
                L = float(line_input[1])
            elif line_input[0] == 'n_cells':
                n_cells = int(line_input[1])
                S_w = zeros(n_cells,float)
                S_n = zeros(n_cells,float)                
                z = zeros(n_cells,float)
                poros = zeros(n_cells,float)
                k = zeros(n_cells,float)
                alpha = zeros(n_cells,float)
                n = zeros(n_cells,float)
                rs = zeros((n_cells,2),float)
                dS_0 = zeros(n_cells,float)
                Ss_0 = zeros(n_cells,float)
                dS_f = zeros(n_cells,float)
                Ss_f = zeros(n_cells,float)
                vol = zeros(n_cells,float)                
            else:    # cross-sectional area of vertical column
                A = float(line_input[1])
        # set up cell objects
        if flag_read[2] and i > 0 and i <= n_cells:                           # note: 1st column contains cell index numbers, for visual reference
            S_w[i-1] = float(line_input[1])
            S_n[i-1] = float(line_input[2])                                     # initial conditions
            z[i-1] = float(line_input[3])
            poros[i-1] = float(line_input[4])
            k[i-1] = float(line_input[5])
            alpha[i-1] = float(line_input[6])
            n[i-1] = float(line_input[7])
            rs[i-1] = [float(line_input[8]),float(line_input[9])]
            dS_0[i-1] = float(line_input[10])
            Ss_0[i-1] = float(line_input[11])
            dS_f[i-1] = float(line_input[12])
            Ss_f[i-1] = float(line_input[13])
            vol[i-1] = float(line_input[14])            
        i += 1
    liquids = Liquid(name,rho,u,beta)
    cells = Cell(z,poros,k,alpha,n,rs,dS_0,dS_f,Ss_0,Ss_f,vol)
    dz = L/n_cells
    input_file.close()
    return liquids,cells,S_w,S_n,n_cells,dz,A

def ReadTimes():
    # read in time steps for which to write output
    t = []
    i = 0
    times_file = open('times.txt','r')
    for line in times_file:
        line_input = line.split()
        if i:t.append(float(line_input[0]))
        i += 1
    times_file.close()
    return array(t)

def WriteOutput(i_step,y_t,t,cells,liquids,n_cells,A,dz,prefix,i_finish):
    # write header
    output_file = open(prefix + '_out_' + str(t[i_step]) + '.txt','w')
    line_out = ['z','\t','P','\t','S_w','\t','S_n','\t','Q_w','\t','Q_n','\n']
    output_file.writelines(line_out)
    # write model results for time t[i_step]
    t_slice = y_t[i_step]
    S_w = t_slice[:n_cells]
    S_n = t_slice[n_cells:]
    for i in xrange(n_cells):
        line_out = [str(cells.z[i]),'\t',str(cells.Pressure(S_w,S_n,liquids)[i]),'\t',str(S_w[i]),'\t',str(S_n[i]),
            '\t',str(cells.Q(0,S_w,S_n,A,dz,liquids)[:-1][i]),'\t',str(cells.Q(1,S_w,S_n,A,dz,liquids)[:-1][i]),'\n']
        output_file.writelines(line_out)
    output_file.close()
    # display saturation results on graph at end-of-run
    if i_finish:
        plt.plot(S_w,cells.z,'-',label = 'Water')
        plt.plot(S_n,cells.z,'-',label = 'NAPL')
        plt.plot(S_w + S_n,cells.z,'-',label = 'Total liquid')        
        plt.xlabel('Saturation')
        plt.ylabel('z')
        plt.legend(loc='upper center',bbox_to_anchor=(0.5, 1.1),ncol=3,fancybox=True,shadow=True)
        plt.show()


###########################################################
#
# main script
#
###########################################################


def VadoseMulti():

    # process input
    prefix = raw_input('Data set prefix for input and output files >')
    liquids,cells,S_w,S_n,n_cells,dz,A = ProcessInput(prefix)
    print 'Read and processed model input file.'

    # time stepping
    t = ReadTimes()
    print 'Read output times.'

    # solve PDEs by method-of-lines
    y_t = cells.ODEHandler(t,S_w,S_n,liquids,n_cells,A,dz)

    # write output
    i_finish = 0
    for i_step in xrange(len(t)):
        if i_step == len(t)-1:i_finish = 1
        WriteOutput(i_step,y_t,t,cells,liquids,n_cells,A,dz,prefix,i_finish)

    print 'Done.'


###################### run script ...

VadoseMulti()  
