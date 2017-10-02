########################################################################
#   Siyanda Matika Absorption Function                                 #
#   18 March 2016                                                      #
#                                                                      #
#                                                                      #
########################################################################
import numpy as np
#import matplotlib.pyplot as plt
#%matplotlib inline
from scipy import integrate
import glob
from utils import find_nearest
from utils import medianArray
import utils as ut




def inter(f,r1,r2):
        s=np.interp(f,r1,r2 )
        return s
def chi(z):
    H0=68.
    oh_m=0.27
    oh_r=9.2e-5
    c=3.0e5
    a =np.sqrt(oh_r*(1.+z)**4. +oh_m*(1.+z)**3.+(1.-oh_m-oh_r))
    b= (H0)*a
    return (c/b)

def n_HI(z,xi,di,i):
    Omega_b_h_2 =0.02230
    var1 = 0.18*(Omega_b_h_2/0.02)*(1.0-xi[i])*(1.0+z)**3.0 * (1.0+di[i] )
    return var1

def Te_(z,xi,di,i):
    return  n_HI(z,xi,di,i)*1.82e-28*(chi(z)*1.0e6*3.08e16) # 1e6*3.08e16 is to convert chi [Mpc] >> m

def A(z_g,z,xi,di,T,i):
    if z>z_g:
        return (1.0 - Te_(z_g,xi,di,i)/T[i] )
    else :
        return 1.0


def r(z):
    return integrate.quad(chi,0,z)[0]


from scipy import integrate as integrate


def A_1(z_g,z,B,K):
    if z_g>z:
        return  1.0-np.exp(-K)*(1.0 + z)**(B)
    else :
        return 1.0



def A_2(z_g,z,a,b,c):
    y= (1.0 + z_g)**(b)*np.exp(a*( np.log(1.0+z_g) )**2.0 +c )#+c
    if z>z_g:
        return ( 1.0-y)
    else :
        return 1.0

def A_2M(z_g,z,a,b,c):
    y= (1.0 + z_g)**(b)*np.exp(a*( np.log(1.0+z_g) )**2.0  )#+c
    if z>z_g:
        return ( 1.0-y)
    else :
        return 1.0


def Snu_Model_A_1(a,B,K,S0,x,z_g):
    nu01 =x*1e-3
    z=(1420.4/x)-1.0
    var1 =1.0#0.214*L0*(46.77)/((1.0+z_g)**(1.0+A) )
    var2 =1.0#(9142.7/r(z_g))**2.0#(1.0e-22/ ( r(z_g)*3.08) )**2.0*(1.0e26*1.0e26*1.0e3/4.0*np.pi)#(9142.7/r(z_g))**2.0
    var3 = S0*(1.4/nu01)**a
    mod  =var1* var2*var3
    model =A_1(z_g,z,B,K)

    return model




def Snu_Model_A(x,B,K,z_g,a,L0):
    """
    observer frame
    xrest/xobs=1+z
    x=freq_obs
    z_g=z_source
    freq_break_
    """
    nu01 =x*1.0e-3
    var1=0.214*(L0/1.0)*(46.77)/( (1.0+z_g)**(1.+a) )
    var2 =(9142.7/ Rnter(z_g) )**2.0#(1.0e-22/ ( r(z_g)*3.08) )**2.0*(1.0/4.0*np.pi)*1.0e26*1.0e26*1.0e3
    var3 = var1*var2*(1.4/nu01)**a#var1*var2*(1.4/nu01)**a
    #xg=1420.4/(1.0+z_g)
    zo=(1420.4/x)-1.0
    if zo<z_g:
        var=  B*np.log(1.0+zo) + K
        return  (1.0-np.exp(var) )*var3
    else:
        return 1.0*var3



def Int_S(x,x_i,B,K,z_g,a,L0):
    var=integrate.quad(lambda nu:(1420.4/nu*nu)*Snu_Model_A(nu,B,K,z_g,a,L0),x,x_i)
    return var[0]#+var[1]


def Snu_Model_A_2(a,A,B,C,S0,x,z_g):
    nu01 =x*1e-3
    z=(1420.4/x)-1.0
    var1 =1.0#0.214*L0*(46.77)/( (1.0+z_g)**(1.0+a) )
    var2 =1.0#(9142.7/r(z_g))**2.0#(1.0e-22/ ( r(z_g)*3.08) )**2.0*(1.0/4.0*np.pi)*1.0e26*1.0e26*1.0e3# (9142.7/r(z_g))**2.0
    var3 = S0*(1.4/nu01)**a
    mod  =var1*var2*var3
    model = mod*A_2(z,z_g,A,B,0.0)
    #A_z(15.0,i,8.0,25.)
    return model

def Snu_Model_A_2_M(a,S0,A,B,x,z_g):
    #L0=1.0
    nu01 =x*1e-3
    z=(1420.4/x)-1.0
    var1 =1.0# 0.214*L0*(46.77)/( (1.0+z_g)**(1.0+a) )
    var2 =1.0#(9142.7/r(z_g))**2.0#(1.0e-22/ ( r(z_g)*3.08) )**2.0*(1.0/4.0*np.pi)*1.0e26*1.0e26*1.0e3# (9142.7/r(z_g))**2.0
    var3 =S0*(1.4/nu01)**a
    mod  =var1*var2*var3
    model = mod*A_2M(z,z_g,A,B,0.0)
    #A_z(15.0,i,8.0,25.)
    return model


def Snu_Model_NoN(a,S0,x):
    nu01 =x*1e-3
    var1 =1.0#0.214*L0*(46.77)/( (1.0+z_g)**(1.0+a) )
    var2 =1.0#(9142.7/r(z_g))**2.0#(1.0e-22/ ( r(z_g)*3.08) )**2.0*(1.0e26*1.0e26*1.0e3/(4.0*np.pi) ) # (9142.7/r(z_g))**2.0
    var3 =S0*(1.4/nu01)**a
    mod  =var1*var2*var3
    model = mod
    return model



def S_n1(x,z_g,L0,a):
    nu01 =x*1.0e-3
    var1=0.214*(L0/1.0)*(46.77)/( (1.0+z_g)**(1.+a) )
    var2 =(9142.7/ r(z_g))**2.0#(1.0e-22/ ( r(z_g)*3.08) )**2.0*(1.0/4.0*np.pi)*1.0e26*1.0e26*1.0e3
    var3 = (1.4/nu01)**(a)
    var4 = 1.0
    return var1* var2*var3*var4


def S_Eor(x,z_g,L0,a,XI,di,T,i):
    z=(1420.4/x)-1.0
    return A(z,z_g,XI,di,T,i)*S_n1(x,z_g,L0,a)


def Int_S_EoR(x,x_i,a,z_g,L0,XI,di,T,i):
    b=integrate.quad(lambda nu:S_Eor(nu,z_g,L0,a,XI,di,T,i),x,x_i)[0]
    return b/(x_i-x)

def Int_S_N(x,x_i,a,z_g,L0):
    b=integrate.quad(lambda nu:S_n1(nu,z_g,L0,a),x,x_i)[0]
    return b



def Tsys(nu):
    '''
    Tsys: Represents the sky temperature at a given frequency
    T_sky: Sky temperature :units [k ]

    args:
    The only argument required is frequncy (nu)

    '''
    T_sky = 60.0*(300./nu)**(2.55)
    T=(0.1*T_sky+40.0)+( T_sky)
    return T

def sig_s(nu,dnu,t_int):
    var1= Tsys(nu)*2.0*1.38e-23
    if nu>110.:
        A_tt = 0.5*866.0*(925.0)*(110./nu)**2.0
    else:
        A_tt =0.5*866.0*(925.0)
    #A_tt = 0.5*866.0*(925.0)*(110./nu)**2.0
    var2= var1*1.0e26/ ( A_tt*np.sqrt(dnu*1.0e6*t_int*60.0*60.0) )
    return var2*1.0e3 #!!NB the 10 dividing is to improve the signal to noise











class FLUX:

	"""
	Master Function which computes the Flux spectrum of qasars in the
	Epoch of Reionization. The data used to generate the flux is
	from SimFast21 simulation by Santos et al (2008).






	"""
	RR=np.loadtxt("Viro7_r_z6_27_25_Ie4.dat",skiprows=1,unpack=True,dtype=float)
	z_eor      =    RR[0,:][~np.isnan( RR[0,:])]

	xi_eor     =     RR[1,:][~np.isnan( RR[0,:])]
	xi1_eor    =     RR[2,:][~np.isnan( RR[0,:])]
	xi2_eor    =     RR[3,:][~np.isnan( RR[0,:])]
	xi3_eor    =     RR[4,:][~np.isnan( RR[0,:])]
	xi4_eor    =     RR[5,:][~np.isnan( RR[0,:])]

	delta_eor  =      RR[6,:][~np.isnan(  RR[0,:])]
	delta1_eor =      RR[7,:][~np.isnan(  RR[0,:])]
	delta2_eor =      RR[8,:][~np.isnan(  RR[0,:])]
	delta3_eor =      RR[9,:][~np.isnan(  RR[0,:])]
	delta4_eor =      RR[10,:][~np.isnan( RR[0,:])]

	Ts_eor          =    RR[11,:][~np.isnan( RR[0,:])]
	Ts1_eor         =    RR[12,:][~np.isnan( RR[0,:])]
	Ts2_eor         =    RR[13,:][~np.isnan( RR[0,:])]
	Ts3_eor         =    RR[14,:][~np.isnan( RR[0,:])]
	Ts4_eor         =    RR[15,:][~np.isnan( RR[0,:])]




	"""
	On  __init__ We initialise our Master function for flux during EoR
	"""
	def __init__(self,red_g,L_g,index_g,nu_noise_space=None,plotting=None):
		self.red_g=red_g;self.L_g=L_g;self.index_g=index_g;

		if nu_noise_space is None:
			self.nu_i=50.0;self.nu_f=250.0;self.N_nu=50.0;self.t_int=1000.0
		else:
			self.nu_i=nu_noise_space["nu_i"]
			self.nu_f=nu_noise_space["nu_f"]
			self.N_nu=nu_noise_space["N_nu"]
			self.t_int=nu_noise_space["t_int"]

		self.nu_arr=self.nu_i+(self.nu_f-self.nu_i)*np.arange(self.N_nu)/(self.N_nu-1.0)
		self.z_arr= (1420.4/self.nu_arr)-1.0
		self.Dnu=(self.nu_f-self.nu_i)/self.N_nu
		#xi
		self.xi_z = inter(self.z_arr,self.z_eor,self.xi_eor)
		self.xi_z1 = inter(self.z_arr,self.z_eor,self.xi1_eor)
		self.xi_z2 = inter(self.z_arr,self.z_eor,self.xi2_eor)
		self.xi_z3 = inter(self.z_arr,self.z_eor,self.xi3_eor)
		self.xi_z4 = inter(self.z_arr,self.z_eor,self.xi4_eor)
		#d_z
		self.d_z = inter(self.z_arr,self.z_eor,self.delta_eor)
		self.d_z1 = inter(self.z_arr,self.z_eor,self.delta1_eor)
		self.d_z2 = inter(self.z_arr,self.z_eor,self.delta2_eor)
		self.d_z3 = inter(self.z_arr,self.z_eor,self.delta3_eor)
		self.d_z4 = inter(self.z_arr,self.z_eor,self.delta4_eor)
		#Ts
		self.Ts_z  = inter(self.z_arr,self.z_eor,self.Ts_eor)
		self.Ts_z1 = inter(self.z_arr,self.z_eor,self.Ts1_eor)
		self.Ts_z2 = inter(self.z_arr,self.z_eor,self.Ts2_eor)
		self.Ts_z3 = inter(self.z_arr,self.z_eor,self.Ts3_eor)
		self.Ts_z4 = inter(self.z_arr,self.z_eor,self.Ts4_eor)

	def formula(self):
		return "v0*t - 0.5*g*t**2; v0=%g" %len(self.xi_z)

	def S_P(self,z,L0,a):
                xm=medianArray(self.nu_arr)#defualt is los 4: 1st try 2,3,noda()
                Sn011  = np.array( [ Int_S_EoR(xm[i],xm[i+1],a,z,L0,self.xi_z,self.d_z,self.Ts_z,i)  for i in range(len(xm[:-1])) ]  )

		Sn1_av=(Sn011)
		return Sn1_av

	def S_non(self):
           xm=medianArray(self.nu_arr)
           z,L0,a=self.red_g,self.L_g,self.index_g
           y=np.ones(len(xm)-1)
           for i,xi in enumerate(xm[:-1]):
              #y[i]=Int_S(self.nu_arr[i],self.nu_arr[i+1],8.0,23.0,z_g)
              y[i]=Int_S_N(xm[i],xm[i+1],a,z,L0)

           return y

	def S_N(self,z,L0,a,Dnu,t_int):
                xm=medianArray(self.nu_arr)
                ssg= np.array([ sig_s(jj,Dnu,t_int) for jj  in xm[:-1] ])
                #tt   = np.array( [ Int_ST(xm[i],xm[i+1],a,z,L0,self.xi_z2,self.d_z2,self.Ts_z2,i)  for i in range(len(xm[:-1])) ]  )

                return (self.nu_arr),self.S_P(z,L0,a),ssg#,tt

	def Result(self):
		z,L0,a=self.red_g,self.L_g,self.index_g
		t,Dnu=self.t_int,self.Dnu

		return self.S_N(z,L0,a,Dnu,t)

	def Con(self):
                z,L0,a=self.red_g,self.L_g,self.index_g
		t,Dnu=self.t_int,self.Dnu
                xm=medianArray(self.nu_arr)
                ssg= np.array([ sig_s(jj,Dnu,t) for jj  in xm[:-1] ])
                Sn014  = np.array( [ Int_S_EoR(xm[i],xm[i+1],a,z,L0,self.xi_z4,self.d_z4,self.Ts_z4,i)  for i in range(len(xm[:-1])) ]  )
                #Sn014  =  np.array( [ S_n1(xm[i],z,L0,a)*A(z,(1420.4/xm[i])-1.0,self.xi_z4,self.d_z4,self.Ts_z4,i) for i in range(len(xm[:-1])) ]  )
                return (self.nu_arr),Sn014,ssg



	def Mod(self,A,B,z,L0,a):
		t,Dnu=self.t_int,self.Dnu
                xm=medianArray(self.nu_arr)
                ssg= np.array([ sig_s(jj,Dnu,t) for jj  in xm[:-1] ])
                #Sn014  =  np.array( [ S_n1(xm[i],z,L0,a)*Snu_Model_A(xm[i],A,B,z) for i in range(len(xm[:-1])) ]  )
                #Sn014  =  np.array( [ Int_S(xm[i],xm[i+1],A,B,z,a,L0)   for i in range(len(xm[:-1])) ]  )
                Sn014  =  np.array( [ Int_S(xm[i],xm[i+1],A,B,z,a,L0)   for i in range(len(xm[:-1])) ]  )

                return Sn014




#nu_space={"nu_i":70.,"nu_f":750,"N_nu":2000.,"Dnu":0.01,"t_int":2000.}
#plotting=False


from utils import medianArray


def Snu_Model_A1(x,B,K,z_g):
    zo=(1420.4/x)-1.0
    if zo<z_g:
        var=  B*np.log(1.0+zo) + K
        return  (1.0-np.exp(var) )
    else:
        return 1.0

def Int_S_1(x,x_i,B,K,z_g):
    var=integrate.quad(lambda nu:Snu_Model_A1(nu,B,K,z_g),x,x_i)[0]#integrate.quad(lambda nu:Snu_Model_A(nu,B,K,z_g),x,x_i)[1]
    return var#+var[1]

def Snu_Model_A_2(x,a,b,c,z_g):
    zo=(1420.4/x)-1.0
    if zo<z_g:
        var=  a*( np.log(1.0+zo) )**2.0+b*np.log(1.0+zo)+ c
        return  (1.0-np.exp(var) )
    else:
        return 1.0

def Int_S2(x,x_i,a,b,c,z_g):
    var=integrate.quad(lambda nu:Snu_Model_A_2(nu,a,b,c,z_g),x,x_i)[0]
    return var

#########################################################

#zr=np.linspace(0,20,10000)
#ru=np.array( [ r(i) for  i in zr ] )
#ru=[r(12.0+i/10.) for i in range(100)]

#from scipy.interpolate import interp1d
#ff=np.loadtxt("rz.txt",dtype=float)
#x1=ff[:,0]
#y1=ff[:,1]

#f = interp1d(x1, y1)
#ru=np.array( [ r(i) for  i in zr ] )

#np.savetxt('rz.txt',np.transpose([zr,ru]))
#plt.plot(x1,y1,label="original")
#plt.plot(zr,ru,"r.",label="uu")
#plt.legend()
#plt.show()
import time
from scipy import stats
start = time.time()
from collections import OrderedDict

class Memoize:
    def __init__(self, f):
        self.f = f
        self.memo = {}
        for x in range(141):
            self.memo[6.0+x/10.] = self.f(6.0+x/10.)
        self.memo=OrderedDict ( sorted(self.memo.items()) )
    def __call__(self, *args):
        #if not args in self.memo:
            #self.memo[args] = self.f(*args)
            #interpolate
        return self.memo[args]

global bb,cc
bb,cc = zip(* Memoize(r).memo.items())

def Rnter(f):
        s=np.interp(f,bb,cc)
        return s
