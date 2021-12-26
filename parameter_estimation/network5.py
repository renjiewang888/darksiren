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
import emcee
import corner
import seaborn as sns
#plt.switch_backend('agg')


#read the data                                                     #nine-dimensional matrix(M_c,eta,iota,d_L,alpha,delta,t_c,phi_c,psi)
datanum=4
fishermatrix1=np.loadtxt("%d/Hdata%d_net.txt"%(datanum,datanum))
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

BHnum=1197          #867, 1197
nn=1000
PH0total=np.ones(nn)


begin_time=time()

while(BHnum<1198):
    fishermatrix=fishermatrix1[[BHnum*9,BHnum*9+1,BHnum*9+2,BHnum*9+3,BHnum*9+4,BHnum*9+5,BHnum*9+6,+BHnum*9+7,BHnum*9+8]]    #read the  fishermartix of ith BH events
    #print(fishermatrix)
    print(BHnum)
    #cov2=np.linalg.inv(fishermatrix)
    #print(cov2)
    #print(cov2[2][2],cov2[4][4],cov2[5][5],cov2[8][8])
    #print(np.sqrt(cov2[3][3]))    
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
    logd_l=math.log(GWd_l,10)

    #separate the position part from the whole covariance matrix
    cov_3para=cov[[3,4,5]]
    cov_3para=cov_3para[:,[3,4,5]]
    #print(cov_3para)
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
    print(alphamin,alphamax,deltamin,deltamax)
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
    if(N>10**7):
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
    #sample galaxies and locate the host galaxies
    with open("%d/galaxy/galaxy%d_2.txt"%(datanum,BHnum),"r+") as fa:
        fa.truncate()
    
    with open("%d/source/source%d_2.txt"%(datanum,BHnum),"r+") as s:
        s.truncate()
    
    
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
        
        
        with open("%d/galaxy/galaxy%d_3.txt"%(datanum,BHnum),"a") as fa:
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
            
            with open("%d/source/source%d_3.txt"%(datanum,BHnum),"a") as s:
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
    
    with open("%d/source/source%d_3.txt"%(datanum,BHnum),"a") as s:
               s.write('%f '%alpha)
               s.write('%f '%delta)
               s.write('%f '%redshift)
    '''
    
    alpha_g,delta_g,d_lg=np.loadtxt("%d/galaxy/galaxy%d.txt"%(datanum,BHnum),usecols=(0,1,3),unpack=True)
    source_alpha,source_delta,source_zz=np.loadtxt("%d/source/source%d.txt"%(datanum,BHnum),usecols=(0,1,2),unpack=True)
    source_logdl=np.zeros(len(source_zz))
    for i in range(len(source_zz)):
        source_logdl[i]=math.log(D_l(source_zz[i])/1000,10)
        source_logdl[i]=np.around(source_logdl[i],decimals=6,out=None)    
    #print(alpha_g.shape,delta_g.shape,d_lg.shape)
    para=np.vstack((d_lg,alpha_g,delta_g))
    source_para=np.vstack((source_logdl,source_alpha,source_delta))
    #print(source_para,para.shape)
    xyz=np.dot(np.transpose(eigenvector),para)
    source_xyz=np.dot(np.transpose(eigenvector),source_para)  

    

    source_x=[]
    source_y=[]
    source_z=[]
    '''
    for i in range(len(source_zz)):
        for j in range(i,len(alpha_g)):
            if(xyz[0][j]==source_xyz[0][i]):
                source_x.append(source_xyz[0][i])
                source_y.append(source_xyz[1][i])
                source_z.append(source_xyz[2][i])
            else:
                galaxy_x.append(xyz[0][j])
                galaxy_y.append(xyz[1][j])
                galaxy_z.append(xyz[2][j]) 
    source_x.append(source_xyz[0][len(source_zz)-1])
    source_y.append(source_xyz[1][len(source_zz)-1])
    source_z.append(source_xyz[2][len(source_zz)-1])    
    '''     
    galaxy_x.append(xyz[0])
    galaxy_y.append(xyz[1])
    galaxy_z.append(xyz[2])  
    source_x.append(source_xyz[0])
    source_y.append(source_xyz[1])
    source_z.append(source_xyz[2])
    
    sns.set_style("white")
           
    fig=plt.figure()  
    
    plt.gcf().set_size_inches(8/2.54,6.2/2.54)             # (14，10)  （10，7） 8/2.54,6.2/2.54
    axes3d=Axes3D(fig)

    logd_l=np.around(logd_l,decimals=6,out=None)
    alpha=np.around(alpha,decimals=6,out=None)
    delta=np.around(delta,decimals=6,out=None)
    vertex=np.array([[logdlmin,logdlmin,logdlmin,logdlmin,logdlmax,logdlmax,logdlmax,logdlmax,logd_l],
                    [alphamin,alphamax,alphamax,alphamin,alphamin,alphamax,alphamax,alphamin,alpha],
                    [deltamin,deltamin,deltamax,deltamax,deltamin,deltamin,deltamax,deltamax,delta]])
    #print(vertex)
    vertex_xyz=np.dot(np.transpose(eigenvector),vertex)
    #print(vertex_xyz)
    
    for i in range(4):
        if(i<3):
            if(i==5):    
                axes3d.plot((vertex_xyz[0][i],vertex_xyz[0][i+1]),(vertex_xyz[1][i],vertex_xyz[1][i+1]),(vertex_xyz[2][i],vertex_xyz[2][i+1]),color='black',linestyle='--',linewidth=0.35)
            else:
                axes3d.plot((vertex_xyz[0][i],vertex_xyz[0][i+1]),(vertex_xyz[1][i],vertex_xyz[1][i+1]),(vertex_xyz[2][i],vertex_xyz[2][i+1]),color='black',linestyle='--',linewidth=0.35)
        else:
            axes3d.plot((vertex_xyz[0][i],vertex_xyz[0][0]),(vertex_xyz[1][i],vertex_xyz[1][0]),(vertex_xyz[2][i],vertex_xyz[2][0]),color='black',linestyle='--',linewidth=0.35)

    for i in range(4,8):
        if(i<7):
            if(i==0):    #6  0
                axes3d.plot((vertex_xyz[0][i],vertex_xyz[0][i+1]),(vertex_xyz[1][i],vertex_xyz[1][i+1]),(vertex_xyz[2][i],vertex_xyz[2][i+1]),color='black',linestyle='--',linewidth=0.35)
            else:
                axes3d.plot((vertex_xyz[0][i],vertex_xyz[0][i+1]),(vertex_xyz[1][i],vertex_xyz[1][i+1]),(vertex_xyz[2][i],vertex_xyz[2][i+1]),color='black',linestyle='--',linewidth=0.35)
        else:
            axes3d.plot((vertex_xyz[0][i],vertex_xyz[0][4]),(vertex_xyz[1][i],vertex_xyz[1][4]),(vertex_xyz[2][i],vertex_xyz[2][4]),color='black',linestyle='--',linewidth=0.35)

    for i in range(4):
        if(i==2):
            axes3d.plot((vertex_xyz[0][i],vertex_xyz[0][i+4]),(vertex_xyz[1][i],vertex_xyz[1][i+4]),(vertex_xyz[2][i],vertex_xyz[2][i+4]),color='black',linestyle='--',linewidth=0.35)
        else:
            axes3d.plot((vertex_xyz[0][i],vertex_xyz[0][i+4]),(vertex_xyz[1][i],vertex_xyz[1][i+4]),(vertex_xyz[2][i],vertex_xyz[2][i+4]),color='black',linestyle='--',linewidth=0.35)
    
    #axes3d.scatter(vertex_xyz[0],vertex_xyz[1],vertex_xyz[2])

    
    #axes3d.scatter(source_x,source_y,source_z,s=20,c='r',marker='*',alpha=1)
    axes3d.scatter(galaxy_x,galaxy_y,galaxy_z,marker=',',s=0.1,c='r',alpha=0.3)        #  5,0.6
    
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
    #axes3d.plot_surface(x,y,z,color='r',rstride=10,cstride=10,linewidth=0.3,edgecolors='black',alpha=0.1)
    axes3d.plot_wireframe(x,y,z,rstride=5, cstride=10,alpha=1,linewidth=0.35)             #0.5
    axes3d.view_init(elev=10,azim=67)        # 10,70  40,65

    plt.grid(linewidth=0.35)
    #plt.legend()
    font={'family':'Arial','weight':'normal','size':5.5}
    axes3d.set_xlabel('x',font,labelpad=0.0)
    axes3d.set_ylabel('y',font,labelpad=0.0)
    axes3d.set_zlabel('z',font,labelpad=0.0)
    plt.tick_params(labelsize=5.5,pad=0.005)          #     14,  12
    #plt.title("galaxy localization")
    axes3d.spines['bottom'].set_linewidth(0.35)
    axes3d.spines['left'].set_linewidth(0.35)
    #axes3d.spines['right'].set_linewidth(0.35)
    axes3d.set_zlim(np.min(galaxy_z)-0.001,np.max(galaxy_z)+0.001)
    axes3d.set_xlim(np.min(galaxy_x)-0.0015,np.max(galaxy_x)+0.0015)
    axes3d.set_ylim(np.min(galaxy_y)-0.0012,np.max(galaxy_y)+0.0012)
    plt.savefig("%d/ellipsoid1/ellipsoid%d_6.pdf"%(datanum,BHnum),bbox_inches='tight',pad_inches=0.0)
    #plt.clf()
    plt.show()
    
    '''
    fig=plt.figure()
    plt.gcf().set_size_inches(15, 10)  
    axes3d=Axes3D(fig)
    aa=eigenvalue[0]*chi_square
    bb=eigenvalue[1]*chi_square
    cc=eigenvalue[2]*chi_square       #long half axis, mid half axis and short half axis of ellipsoid
    u=np.linspace(0,2*np.pi,100)    
    v=np.linspace(0,np.pi,100)
    A=mean_xyz[0]+np.sqrt(aa)*np.outer(np.cos(u),np.sin(v))
    B=mean_xyz[1]+np.sqrt(bb)*np.outer(np.sin(u),np.sin(v))
    C=mean_xyz[2]+np.sqrt(cc)*np.outer(np.ones(np.size(u)),np.cos(v))     
    #print(A[0,0])
    i=0
    j=0
    k=0
    x=np.ones([100,100])
    y=np.ones([100,100])
    z=np.ones([100,100])
    points=np.ones([10000,3])
    while(i<100):
        while(j<100):
            coefficient=np.array([[A[i,j]],[B[i,j]],[C[i,j]]])
            #print(coefficient)
            #print(np.transpose(eigenvector))
            w=np.linalg.solve(np.transpose(eigenvector),coefficient)          #(logdl,alpha,delta)
            #print(w)                                 
            points[k][0]=x[i][j]=10**(w[0])*np.cos(w[2])*np.sin(w[1])*1000
            points[k][1]=y[i][j]=10**(w[0])*np.sin(w[2])*np.sin(w[1])*1000
            points[k][2]=z[i][j]=10**(w[0])*np.cos(w[1])*1000
            k=k+1
            j=j+1
        j=0
        i=i+1
    axes3d.plot_surface(x,y,z,color='r')
    axes3d.view_init(elev=30,azim=65)
    axes3d.set_xlabel('x/Mpc')
    axes3d.set_ylabel('y/Mpc')
    axes3d.set_zlabel('z/Mpc')
    #axes3d.set_xlim(np.min(x)-500,0)
    #axes3d.set_ylim(0,np.max(y)+500)
    #axes3d.set_zlim(np.min(z),0)
    plt.savefig("%d/ellipsoid1/ellipsoid_%d"%(datanum,BHnum),dpi=300)
    #plt.clf()
    #plt.show()
    
    
    
    
    from matplotlib import animation 
    def rotate(angle): 
        axes3d.view_init(azim=angle)
    rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0,362,4),interval=100)
    rot_animation.save('rotation.gif', dpi=80)
    '''
    
    if(N>0):
        source_alpha,source_delta,source_zz=np.loadtxt("%d/source/source%d.txt"%(datanum,BHnum),usecols=(0,1,2),unpack=True)
        num=len(source_alpha)-1
    else:
        source_a,source_d,source_z=np.loadtxt("%d/source/source%d.txt"%(datanum,BHnum),usecols=(0,1,2),unpack=True)
        source_alpha=[source_a]
        source_delta=[source_d]
        source_zz=[source_z]
        num=0
    print(num)
    '''
    with open("%d/result%d.txt"%(datanum,datanum),"a") as fa:
               fa.write('%d\n'%BHnum)
               #fa.write('%f '%redshift)
               #fa.write('%f '%M)
               #fa.write('%f '%iota)
               #fa.write('%f '%alpha)
               #fa.write('%f '%delta) 
               #fa.write('%f '%omega)
               #fa.write('%f '%cov_3para[0][0])
               #fa.write('%d \n'%num)
    
    '''
    
    
    if(num<50):
        nn=10000
        h=np.linspace(65,71,nn)
    elif(num<50000):
          nn=1000
          h=np.linspace(40,100,nn)
    else:
         nn=100
         h=np.linspace(30,120,nn)

    
    
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
    
    '''
    if(((h[uppervalue]-h[inx])<67.74*0.05)&((h[uppervalue]-h[inx])>67.74*0.01)):
        for i in range(len(PH0total)):
            PH0total[i]=PH0total[i]*P_H0GW[i]
        plt.gcf().set_size_inches(10, 7)    
        plt.plot(h,P_H0GW,color='grey')  
        plt.xlabel('$H_0$($km\ s^{-1}$$Mpc^{-1}$)')
        plt.ylabel('p($H_0$)($km^{-1}s Mpc$)')
    '''
    plt.plot(h,P_H0GW,label="%.2f^{%.2f}_{%.2f}"%(h[inx],(h[uppervalue]-h[inx]),(h[lowvalue]-h[inx])))
    plt.title("n=%f N=%d"%(n,num))
    plt.legend()
    #plt.ylim((0,P_H0GW[inx]+0.001))
    plt.xlabel('$H_0$[$km\ s^{-1}$$Mpc^{-1}$]')
    plt.ylabel('p($H_0$)')
    #plt.savefig("%d/image1/test/image%d_3"%(datanum,BHnum))
    #plt.clf()
    #plt.close('all')
    plt.show()
    h=np.mat(h)
    P_H0GW=np.mat(P_H0GW)
    #print(h.shape,P_H0GW.shape)
    h_P=np.vstack((h,P_H0GW))
    h_P=np.transpose(h_P)
    #print(h_P.shape)
    #np.savetxt("%d/image1/h_P%dnet.txt"%(datanum,BHnum),h_P)
    

    BHnum=BHnum+1    
    


end_time=time()
run_time=end_time-begin_time
print(run_time)
'''
BHnum=992
h,P_H0GW=np.loadtxt("%d/image1/h_P%dnet.txt"%(datanum,BHnum),usecols=(0,1),unpack=True)
h=np.array(h)
inx=np.argmax(P_H0GW)
print(BHnum,h[inx])
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
for i in range(len(PH0total)):
    PH0total[i]=PH0total[i]*P_H0GW[i]
plt.gcf().set_size_inches(10, 7)    
plt.plot(h,P_H0GW,color='grey')  
plt.xlabel('$H_0$($km\ s^{-1}$$Mpc^{-1}$)')
plt.ylabel('p($H_0$)($km^{-1}s Mpc$)')
BHnum=971
h,P_H0GW=np.loadtxt("%d/image1/h_P%dnet.txt"%(datanum,BHnum),usecols=(0,1),unpack=True)
h=np.array(h)
inx=np.argmax(P_H0GW)
print(BHnum,h[inx])
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
for i in range(len(PH0total)):
    PH0total[i]=PH0total[i]*P_H0GW[i]
plt.plot(h,P_H0GW,color='grey')  
BHnum=1445
h,P_H0GW=np.loadtxt("%d/image1/h_P%dnet.txt"%(datanum,BHnum),usecols=(0,1),unpack=True)
h=np.array(h)
inx=np.argmax(P_H0GW)
print(BHnum,h[inx])
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
for i in range(len(PH0total)):
    PH0total[i]=PH0total[i]*P_H0GW[i]
plt.plot(h,P_H0GW,color='grey')  
BHnum=1203
h,P_H0GW=np.loadtxt("%d/image1/h_P%dnet.txt"%(datanum,BHnum),usecols=(0,1),unpack=True)
h=np.array(h)
inx=np.argmax(P_H0GW)
print(BHnum,h[inx])
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
for i in range(len(PH0total)):
    PH0total[i]=PH0total[i]*P_H0GW[i]
plt.plot(h,P_H0GW,color='grey')  
'''

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

plt.savefig("Q3dnet15%d"%datanum)
#plt.show()
'''

