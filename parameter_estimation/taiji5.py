import numpy as np
import matplotlib
import math
from scipy import integrate
from scipy.optimize import fsolve,root
import sys
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import norm
import sympy
from time import * 
from scipy.spatial import ConvexHull
from tqdm import tqdm
from scipy import interpolate
#plt.switch_backend('agg')


#read the data                                                     #nine-dimensional matrix(M_c,eta,iota,d_L,alpha,delta,t_c,phi_c,psi)
datanum=4
fishermatrix1=np.loadtxt("%d/Hdata%d_taiji.txt"%(datanum,datanum))
#print(fishermatrix1.shape)

#print(fishermatrix)
redshift1=np.loadtxt("%d/redshift%d.txt"%(datanum,datanum))
iota1=np.loadtxt("%d/iota.txt"%datanum)
alpha1=np.loadtxt("%d/theta.txt"%datanum)
delta1=np.loadtxt("%d/phi.txt"%datanum)
psi1=np.loadtxt("%d/psi.txt"%datanum)
M1=np.loadtxt("%d/Mc%d.txt"%(datanum,datanum))     #Msun



print('maximum redshift',np.max(redshift1))
print('minimum redshift',np.min(redshift1))

#variables
omega_m=0.3
omega_lambda=0.7
c=3.0*10**5
eta=0.25
            #rad
         #mean of alpha

t_c=0
phi_c=2    #rad
chi_square=11.34              #99%confidence 


def E(z):
  return (omega_m*(1+z)**3+omega_lambda)**(0.5)
def X(z):     #comoving distance
  def E1(z):
	  return ((omega_m*(1+z)**3+omega_lambda)**(-0.5))
  return (c/H_0*integrate.quad(E1,0,z)[0])           #Mpc
def H(z):
	return (H_0*E(z))

def D_l(z):
	return((1+z)*X(z))                 #luminosity distance  

redshifttxt,comovingtxt,logd_ltxt=np.loadtxt('redshift_logdl.txt',usecols=(0,1,2),unpack=True)
f= interpolate.interp1d(comovingtxt,redshifttxt,kind=3)


BHnum=859
nn=1000
PH0total=np.ones(nn)


begin_time=time()

while(BHnum<860):
    fishermatrix=fishermatrix1[[BHnum*9,BHnum*9+1,BHnum*9+2,BHnum*9+3,BHnum*9+4,BHnum*9+5,BHnum*9+6,+BHnum*9+7,BHnum*9+8]]    #read the  fishermartix of ith BH events
    #print(fishermatrix)
    print(BHnum)
    cov2=np.linalg.inv(fishermatrix)
    print(cov2[3][3],cov2[4][4],cov2[5][5])
    print(np.sqrt(cov2[3][3]))
    redshift=redshift1[BHnum]
    iota=iota1[BHnum]
    print('iota',iota)
    alpha=alpha1[BHnum]
    delta=delta1[BHnum]
    psi=psi1[BHnum]
    M=M1[BHnum]
    print('redshift=',redshift,'alpha=',alpha,'delta=',delta,'Mc',M)
    H_0=67.74
    GWd_l=D_l(redshift)/1000               #Gpc
    print('GWd_l=',GWd_l)
    print('X=',X(redshift))
    fishermatrix[3,:]=fishermatrix[3,:]*GWd_l
    fishermatrix[:,3]=fishermatrix[:,3]*GWd_l
    cov=np.linalg.inv(fishermatrix)
    #print(cov)
    omega=2*np.pi*np.sin(alpha)*np.sqrt(cov[4][4]*cov[5][5]-cov[4][5]**2)
    omega=omega/(2*np.pi)**2*360**2
    print('omega',omega,'deg^2')

    #print(np.sqrt(cov[0][0]),np.sqrt(cov[2][2]))
    
    #separate the position part from the whole covariance matrix
    cov_3para=cov[[3,4,5]]
    cov_3para=cov_3para[:,[3,4,5]]
    '''
    cov_3para=np.zeros((3,3))       #cov(logd_L.alpha,delta)
    i=3
    j=3
    a=0
    while(i<6):
        b=0
        while(j<6):
            cov_3para[a][b]=cov[i][j]
            j=j+1
            b=b+1;
        i=i+1
        a=a+1
        j=3 
    '''
    print(cov_3para)
    
    
    eigenvalue,eigenvector=np.linalg.eig(cov_3para)    #eigenvalues and eigenvector of covariance matrix
    print('eigenvalue=',eigenvalue)
    #print(eigenvector)
    cov_xyz=np.dot(np.dot(np.linalg.inv(eigenvector),cov_3para),eigenvector)     #diagonalization of covariance matrix =A**(-1)*C*A
    cov_xyz=np.around(cov_xyz,decimals=10,out=None)
    print(cov_xyz)
    def compute_xyz(d_l,alpha,delta):
        para3=np.matrix([[d_l],[alpha],[delta]])
        xyz=np.transpose(eigenvector)*para3
        return xyz;
    logd_l=math.log(GWd_l,10)
    print('GWlogdl',logd_l)
    mean_xyz=compute_xyz(logd_l,alpha,delta)
    #print(mean_xyz)
    
    #mock galaxy catalogue
    n=0.02
    alphamax=alpha+np.sqrt(cov_3para[1][1])*4
    alphamin=alpha-np.sqrt(cov_3para[1][1])*4
    deltamax=delta+np.sqrt(cov_3para[2][2])*4
    deltamin=delta-np.sqrt(cov_3para[2][2])*4
    logdlmax=logd_l+np.sqrt(cov_3para[0][0])*4
    logdlmin=logd_l-np.sqrt(cov_3para[0][0])*4

    
    print(alphamin,alphamax,deltamin,deltamax,logdlmin,logdlmax)
    if((alphamin<(-np.pi))or(alphamax>(2*np.pi))or((alphamax-alphamin)>=(2*np.pi))):
        alphamin=0
        alphamax=np.pi
    if(deltamin<(-2*np.pi)or(deltamax>(4*np.pi))):
        deltamin=0
        deltamax=2*np.pi
    
    def logdl_comoving(logdl):
        dl=10**(logdl)*1000              #Mpc
        def ff(z):
            return (D_l(z)-dl)
        reshift=fsolve(ff,[int(redshift)])
        result=X(reshift)
        return (result);
    if(logdlmax>19):
        BHnum=BHnum+1
        continue
    comovmax=logdl_comoving(logdlmax)
    comovmin=logdl_comoving(logdlmin)
    print('X min and max',comovmin,comovmax)
    V_max=(deltamax-deltamin)*(comovmax**3-comovmin**3)*(np.cos(alphamin)-np.cos(alphamax))/3
    print('V_max=',V_max,'Mpc^3')
    if(V_max>0):
        N=int(n*V_max)
    else:
        print('false')
        BHnum=BHnum+1
        continue
    print('N=',N)

    if(N>10**6):
        BHnum=BHnum+1
        continue

    if(N==1):
        delta_comoving=0
    else:
         delta_comoving=(comovmax-comovmin)/(N-1)
    num=0
    source_zz=[]     #hold the redshift of the possible source
    source_alpha=[]
    source_delta=[]
    galaxy_x=[]
    galaxy_y=[]
    galaxy_z=[]
    
    
    '''
    with open("%d/taijigalaxy/galaxy%d.txt"%(datanum,BHnum),"r+") as fa:
        fa.truncate()
    with open("%d/taijisource/source%d.txt"%(datanum,BHnum),"r+") as s:
        s.truncate()
    '''
    #sample galaxies and locate the host galaxies
    for i in range(0,N):
    
        alpha_g=np.random.uniform(alphamin,alphamax,1)
        delta_g=np.random.uniform(deltamin,deltamax,1)
        if(alpha_g<0):
            alpha_g=-alpha_g
        if(alpha_g>(np.pi)):
            alpha_g=np.pi-(alpha_g-np.pi)
        if(delta_g<0):
            delta_g=2*np.pi+delta_g
        if(delta_g>(2*np.pi)):
            delta_g=delta_g-2*np.pi
        def fun(comoving):
            def f(z):
                return (X(z)-comoving)
            result=fsolve(f,[int(redshift)])
            return (result)    
        comov=np.random.uniform(comovmin,comovmax,1)
        #comov=comovmin+i*delta_comoving                    #galaxies are uniformly distributed in comoving volume
        comov1=math.log(comov,10)
        if((comov1<1.344757)or(comov1>4.045667)):
             redshift_g=fun(comov1)
        else:
             redshift_g=f(comov1)
        #print('comov',comov,X(redshift_g))
        d_lg=(1.0+redshift_g)*comov/1000
        d_lg=math.log(d_lg,10)
        
        
        with open("%d/taijigalaxy/galaxy%d.txt"%(datanum,BHnum),"a") as fa:
               fa.write('%f '%alpha_g)
               fa.write('%f '%delta_g)
               fa.write('%f '%redshift_g)
               fa.write('%f \n'%d_lg)
        
        
        
        
        para=np.vstack((d_lg,alpha_g,delta_g))
        xyz=np.dot(np.transpose(eigenvector),para) 
        galaxy_x.append(xyz[0][0])
        galaxy_y.append(xyz[1][0])
        galaxy_z.append(xyz[2][0])
        
        #locate the host galaxies
        a=(xyz[0][0]-mean_xyz[0])**2/(eigenvalue[0])+(xyz[1][0]-mean_xyz[1])**2/(eigenvalue[1])+(xyz[2][0]-mean_xyz[2])**2/(eigenvalue[2])
        if(a<=chi_square):
            num=num+1
            source_zz.append(redshift_g)
            source_alpha.append(alpha_g)
            source_delta.append(delta_g)
            
            with open("%d/taijisource/source%d.txt"%(datanum,BHnum),"a") as s:
               s.write('%f '%alpha_g)
               s.write('%f '%delta_g)
               s.write('%f \n'%redshift_g)
            
        #sleep(0.1)
    source_zz.append(redshift)
    source_alpha.append(alpha)
    source_delta.append(delta)
    source_zz=np.array(source_zz)
    source_alpha=np.array(source_alpha)
    source_delta=np.array(source_delta)
    print(num)   
    #print(source_zz)
    #print(source_alpha)
    #print(source_delta) 
    
    with open("%d/taijisource/source%d.txt"%(datanum,BHnum),"a") as s:
               s.write('%f '%alpha)
               s.write('%f '%delta)
               s.write('%f '%redshift)
     
        
    ###
    fig=plt.figure()
    axes3d=Axes3D(fig)
    axes3d.scatter(galaxy_x,galaxy_y,galaxy_z,s=1,c='gray',alpha=0.5)     #散点图
    
    #ellipsoid
    chi_square=11.34              #99%confidence 
    aa=eigenvalue[0]*chi_square
    bb=eigenvalue[1]*chi_square
    cc=eigenvalue[2]*chi_square       #long half axis, mid half axis and short half axis of ellipsoid
    u=np.linspace(0,2*np.pi,100)    
    v=np.linspace(0,np.pi,100)
    x=mean_xyz[0]+np.sqrt(aa)*np.outer(np.cos(u),np.sin(v))
    y=mean_xyz[1]+np.sqrt(bb)*np.outer(np.sin(u),np.sin(v))
    z=mean_xyz[2]+np.sqrt(cc)*np.outer(np.ones(np.size(u)),np.cos(v))     #parametric equation of ellipsoid
    axes3d.plot_surface(x,y,z,color='r')
    axes3d.set_xlabel('x')
    axes3d.set_ylabel('y')
    axes3d.set_zlabel('z')
    plt.savefig("%d/ellipsoidtaiji1/ellipsoid%d"%(datanum,BHnum))
    plt.clf()
    plt.close('all')
    #plt.show()
    
    
    
    
    ####################parameter estimation
    '''
    if(N>0):
        source_alpha,source_delta,source_zz=np.loadtxt("%d/taijisource/source%d.txt"%(datanum,BHnum),usecols=(0,1,2),unpack=True)
        num=len(source_alpha)-1
    else:
        source_a,source_d,source_z=np.loadtxt("%d/taijisource/source%d.txt"%(datanum,BHnum),usecols=(0,1,2),unpack=True)
        source_alpha=[source_a]
        source_delta=[source_d]
        source_zz=[source_z]
        num=0
    print(num)  
    '''

             
               
               
              
    if(num<50):
          nn=5000
          h=np.linspace(64,72,nn)
    elif(num<50000):
          nn=2000
          h=np.linspace(60,80,nn)
    else:
         nn=200
         h=np.linspace(20,130,nn)
    
    h=np.array(h)
    
    
    
    
    ##measure the Hubble constant
    #the prior probability of H_0
    P_H0=1.0/(120-40)        #choose[40,120]km s^-1 Mpc^-1 as prior 
    #the group prior p_0(z), assume the groups are uniformly distributed 
    def p_0(z,h):
        H_0=h
        return(X(z)*X(z)/H(z)/V_max)     #the group prior
    

    # the GW likehood
    
    mu_dl=logd_l
    mu_alpha=alpha
    mu_delta=delta
    fishermatrix3para=np.linalg.inv(cov_3para)
    print(fishermatrix3para)
    C=np.sqrt(2*np.pi)**3*np.sqrt(np.linalg.det(cov_3para))
    def p_GW(z,alpha,delta):
        D_L=D_l(z)/1000
        logD_L=math.log(D_L,10)
        delta_dl=logD_L-mu_dl
        delta_alpha=alpha-mu_alpha
        delta_delta=delta-mu_delta
        #print('logdl',delta_dl,'alpha',delta_alpha,'delta',delta_delta)
        parameter=np.array([delta_dl,delta_alpha,delta_delta])
        result1=np.dot(parameter,fishermatrix3para)
        #print(result1)
        result2=np.dot(result1,np.transpose(parameter))
        result4=delta_dl**2/cov_3para[0][0]
        return(np.exp(-(result2)/2)/C)   #GW  likelihood 
    
    #the posterior the H_0
    def posterior(h):
        H_0=h
        P=0
        i=0
        while(i<(num+1)):
            P=P+P_H0*p_GW(source_zz[i],source_alpha[i],source_delta[i])*p_0(source_zz[i],H_0)
            i=i+1;
        del H_0
        return P
    
    del H_0
    if((num+1)>0):
        P_H0GW=np.zeros(nn)
        j=0
        while(j<nn):
            H_0=h[j]
            P_H0GW[j]=posterior(H_0)
            j=j+1;
        #print(P_H0GW)
        normalization=np.sum(P_H0GW)
        P_H0GW=P_H0GW/normalization
        i=0
        

    inx=np.argmax(P_H0GW)
    print(h[inx])
    if((inx>0)&(inx<nn)):
        possible=P_H0GW[inx]
        uppervalue=inx
        lowvalue=inx
        while((possible<0.68/2)&(uppervalue<(nn-1))):
            uppervalue=uppervalue+1
            possible=possible+P_H0GW[uppervalue];
        print(possible,uppervalue)
        possible=P_H0GW[inx]
        while((possible<0.68/2)&(lowvalue>(0))):
            lowvalue=lowvalue-1
            possible=possible+P_H0GW[lowvalue]
        print(possible,lowvalue)
        print((h[lowvalue]-h[inx]),(h[uppervalue]-h[inx]))

    
    
    plt.plot(h,P_H0GW)
    #plt.ylim((0,P_H0GW[inx]+0.001))
    plt.xlabel('$H_0$[$km\ s^{-1}$$Mpc^{-1}$]')
    plt.ylabel('p($H_0$)')
    plt.savefig("%d/taijiimage1/image%d"%(datanum,BHnum))
    plt.clf()
    plt.close('all')
    #plt.show()
    h=np.mat(h)
    P_H0GW=np.mat(P_H0GW)
    #print(h.shape,P_H0GW.shape)
    h_P=np.vstack((h,P_H0GW))
    h_P=np.transpose(h_P)
    #print(h_P.shape)
    np.savetxt("%d/taijiimage1/h_P%dtj.txt"%(datanum,BHnum),h_P)
    
  
    BHnum=BHnum+1;



end_time=time()
run_time=end_time-begin_time
print(run_time)

'''
#np.savetxt("galaxynum.txt",galaxynum)
normalization1=np.sum(PH0total)
PH0total=PH0total/normalization1
#print(np.sum(PH0total))
max_inx=np.argmax(PH0total)
#print(max_inx)
print('H0',h[max_inx],PH0total[max_inx])
if((max_inx>0)&(max_inx<nn)):
    possible=PH0total[max_inx]
    uppervalue=max_inx
    lowvalue=max_inx
    while((possible<0.68/2)&(uppervalue<(nn-1))):
        uppervalue=uppervalue+1
        possible=possible+PH0total[uppervalue];
    print(possible,uppervalue)
    possible=PH0total[max_inx]
    while((possible<0.68/2)&(lowvalue>(0))):
        lowvalue=lowvalue-1
        possible=possible+PH0total[lowvalue]
    print(possible,lowvalue)
    print((h[lowvalue]-h[max_inx]),(h[uppervalue]-h[max_inx]))

#plt.scatter(h,PH0total)
plt.plot(h,PH0total,color='r')
#plt.xlim(60,80)
plt.ylim((0,PH0total[max_inx]+0.001))
plt.xlabel('$H_0$[$km\ s^{-1}$$Mpc^{-1}$]')
plt.ylabel('p($H_0$)')
#plt.savefig("%d/networkresult01"%datanum)
#plt.show()
'''


