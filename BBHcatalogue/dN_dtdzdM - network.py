import numpy as np
import matplotlib
import math
from scipy import integrate
from scipy.optimize import fsolve,root
import sys
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from matplotlib import ticker
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import norm
from matplotlib import colors
np.set_printoptions(suppress=False)
import seaborn as sns

omega_m=0.32
omega_lambda=0.68
c=2.99792458*10**5     #km/s
Mpckm=3.08359979628196608*10**16*10**3
cMpcyr=c/Mpckm*365*24*60*60     #Mpc/yr
print(cMpcyr)
H_0=67.74

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


data=np.fromfile('Klein16_PopIII.dat',dtype=float)



redshift1,m11,m21,W_PS1=np.loadtxt("Klein16_PopIII.dat",usecols=(0,1,2,22),unpack=True)
redshift2,m12,m22,W_PS2=np.loadtxt("Klein16_Q3delays.dat",usecols=(0,1,2,22),unpack=True)
redshift3,m13,m23,W_PS3=np.loadtxt("Klein16_Q3nodelays.dat",usecols=(0,1,2,22),unpack=True)
print(np.max(redshift1),np.max(redshift2),np.max(redshift3))

#print(np.max(redshift),np.min(redshift))
#print(np.max(m1),np.min(m1),np.max(m2),np.min(m2))
#print(W_PS1)
a1=redshift1.shape[0]
a2=redshift2.shape[0]
a3=redshift3.shape[0]
print(a1,a2,a3)

#merger rates per unit redshift
delta_z=1
N=int(20/delta_z)
print(N)

z=np.zeros(N)
for i in range(len(z)):
    z[i]=i*delta_z
  
#np.savetxt("redshift.txt",z,fmt='%d')

#for i in range(0,20):
	#z.append(i+0.5)

dN_dtdz1=np.zeros(N)
dN_dtdz2=np.zeros(N)
dN_dtdz3=np.zeros(N)
i=0
while(i<a1):
    inx=int(redshift1[i]/delta_z)
    dN_dtdz1[inx]=dN_dtdz1[inx]+4*np.pi*cMpcyr*W_PS1[i]*X((redshift1[i]))**2/delta_z
    i=i+1;
i=0
while(i<a2):
    inx=int(redshift2[i]/delta_z)
    dN_dtdz2[inx]=dN_dtdz2[inx]+4*np.pi*cMpcyr*W_PS2[i]*X((redshift2[i]))**2/delta_z
    i=i+1;
i=0
while(i<a3):
    inx=int(redshift3[i]/delta_z)
    dN_dtdz3[inx]=dN_dtdz3[inx]+4*np.pi*cMpcyr*W_PS3[i]*X((redshift3[i]))**2/delta_z
    i=i+1;

print('z<3,popIII',dN_dtdz1[0]+dN_dtdz1[1]+dN_dtdz1[2])
print('z<3,Q3d',dN_dtdz2[0]+dN_dtdz2[1]+dN_dtdz2[2])
print('z<3,Q3nod',dN_dtdz3[0]+dN_dtdz3[1]+dN_dtdz3[2])

'''
plt.bar(z,height=dN_dtdz1,width=0.3,color='orange',alpha=0.5,label="PopIII")
plt.bar(z,height=dN_dtdz2,width=0.3,color='red',label="Q3d")
plt.bar(z,height=dN_dtdz3,width=0.3,color='green',alpha=0.5,label="Q3nod")

plt.xlabel('redshift')
plt.ylabel('dN/dtdz[$yr^{-1}$]')

plt.plot(z,dN_dtdz1,color='red',label='PopIII')
plt.plot(z,dN_dtdz2,color='orange',label='Q3d')
plt.plot(z,dN_dtdz3,color='green',label='Q3nod')

ax1=plt.gca()     #返回坐标轴
ax1.legend(loc='upper right', fontsize=12, frameon=True, fancybox=True, framealpha=0.2, borderpad=0.3,
           ncol=1, markerfirst=True, markerscale=1, numpoints=1, handlelength=3.5)
plt.show()
'''
#i=0
#dN_dtdz1=np.zeros(N)
#for i in range(len(z)):
    #j=0
    #while(j<a):
        #if((redshift[j]>=i*delta_z)&(redshift[j]<(i+1)*delta_z)):
            #dN_dtdz1[i]=dN_dtdz1[i]+4*np.pi*cMpcyr*W_PS[j]*X(redshift[j])**2/delta_z
        #j=j+1;
    #i=i+1;
#print(dN_dtdz1)
#plt.bar(z,height=dN_dtdz1,width=0.3)

#plt.show()
        
#print(dN_dtdz[0],dN_dtdz1[0])


#merger rates per unit total redshifted mass
M_z1=np.zeros(a1)
M_z2=np.zeros(a2)
M_z3=np.zeros(a3)
i=0
while(i<a1):
    M=m11[i]+m21[i]
    M_z1[i]=M*(m11[i]*m21[i]/M**2)**(3/5)
    #M_z1[i]=m11[i]+m21[i]
    M_z1[i]=math.log((M_z1[i]*(1+redshift1[i])),10)
    i=i+1;
i=0
while(i<a2):
    M=m12[i]+m22[i]
    M_z2[i]=M*(m12[i]*m22[i]/M**2)**(3/5)
    #M_z2[i]=m12[i]+m22[i]
    M_z2[i]=math.log((M_z2[i]*(1+redshift2[i])),10)
    i=i+1;
i=0
while(i<a3):
    M=m13[i]+m23[i]
    M_z3[i]=M*(m13[i]*m23[i]/M**2)**(3/5)
    #M_z3[i]=m13[i]+m23[i]
    M_z3[i]=math.log((M_z3[i]*(1+redshift3[i])),10)
    i=i+1;
print(np.min(M_z1),np.min(M_z2),np.min(M_z3),np.max(M_z1),np.max(M_z2),np.max(M_z3))


delta_M=0.45
num=int((11-2)/delta_M)
print(num)

log_M=np.zeros(num)
i=0
for i in range(len(log_M)):
    log_M[i]=2+i*delta_M

print(log_M)
#np.savetxt("logM.txt",log_M,fmt='%f')
   
'''
dN_dtdM1=np.zeros(num)
dN_dtdM2=np.zeros(num)
dN_dtdM3=np.zeros(num)

i=0
while(i<a1):
    inx=int((M_z1[i]-2.0)/delta_M)
    dN_dtdM1[inx]=dN_dtdM1[inx]+4*np.pi*cMpcyr*W_PS1[i]*X(redshift1[i])**2/delta_M
    i=i+1;
i=0
while(i<a2):
    inx=int((M_z2[i]-2.0)/delta_M)
    dN_dtdM2[inx]=dN_dtdM2[inx]+4*np.pi*cMpcyr*W_PS2[i]*X(redshift2[i])**2/delta_M
    i=i+1;
i=0
while(i<a3):
    inx=int((M_z3[i]-2.0)/delta_M)
    dN_dtdM3[inx]=dN_dtdM3[inx]+4*np.pi*cMpcyr*W_PS3[i]*X(redshift3[i])**2/delta_M
    i=i+1;

i=0
while(i<num):
    if(dN_dtdM1[i]!=0):
        dN_dtdM1[i]=math.log(dN_dtdM1[i],10)
    else:
        dN_dtdM1[i]=-1000
    i=i+1;
i=0
while(i<num):
    if(dN_dtdM2[i]!=0):
        dN_dtdM2[i]=math.log(dN_dtdM2[i],10)
    else:
        dN_dtdM2[i]=-1000
    i=i+1;
i=0
while(i<num):
    if(dN_dtdM3[i]!=0):
        dN_dtdM3[i]=math.log(dN_dtdM3[i],10)
    else:
        dN_dtdM3[i]=-1000
    i=i+1;


#plt.bar(log_M,height=dN_dtdM1,width=delta_M,color='orange',alpha=0.5,label="PopIII")
#plt.bar(log_M,height=dN_dtdM2,width=delta_M,color='red',label="Q3delays")
#plt.bar(log_M,height=dN_dtdM3,width=delta_M,color='green',alpha=0.5,label="Q3nodelays")
plt.xlabel('$\log{[M_c/M_{\odot}]}$')
plt.ylabel('dN/dtdlogM_c[$yr^{-1}$]')
plt.xlim((2,11))
plt.ylim((-1,3.0))
plt.yticks([-1,0,1,2],[0.1,1,10,100])
#plt.xticks([2,3,4,5,6,7,8,9,10,11],['$10^2$','$10^3$','$10^4$','$10^5$','$10^6$','$10^7$','$10^8$','$10^9$','$10^{10}$','$10^{11}$'])
plt.plot(log_M,dN_dtdM1,color='red',label='PopIII')
plt.plot(log_M,dN_dtdM2,color='orange',label='Q3d')
plt.plot(log_M,dN_dtdM3,color='green',label='Q3nod')
ax2=plt.gca()     #返回坐标轴
ax2.legend(loc='upper right', fontsize=12, frameon=True, fancybox=True, framealpha=0.2, borderpad=0.3,
           ncol=1, markerfirst=True, markerscale=1, numpoints=1, handlelength=3.5)
plt.show()
'''
#merger rates per unit total redshifted mass and per unit redshift
dN_dtdzdM1=np.zeros((num,N))
dN_dtdzdM2=np.zeros((num,N))
dN_dtdzdM3=np.zeros((num,N))
dN_dt1=np.zeros((num,N))
dN_dt2=np.zeros((num,N))
dN_dt3=np.zeros((num,N))
i=0
while(i<a1):
    inx1=int(redshift1[i]/delta_z)
    inx2=int((M_z1[i]-2)/delta_M)
    dN_dtdzdM1[inx2][inx1]=dN_dtdzdM1[inx2][inx1]+4*np.pi*cMpcyr*W_PS1[i]*X(redshift1[i])**2/delta_M/delta_z
    dN_dt1[inx2][inx1]=dN_dt1[inx2][inx1]+4*np.pi*cMpcyr*W_PS1[i]*X(redshift1[i])**2
    i=i+1;
i=0
while(i<a2):
    inx1=int(redshift2[i]/delta_z)
    inx2=int((M_z2[i]-2)/delta_M)
    dN_dtdzdM2[inx2][inx1]=dN_dtdzdM2[inx2][inx1]+4*np.pi*cMpcyr*W_PS2[i]*X(redshift2[i])**2/delta_M/delta_z
    dN_dt2[inx2][inx1]=dN_dt2[inx2][inx1]+4*np.pi*cMpcyr*W_PS2[i]*X(redshift2[i])**2
    i=i+1;
i=0
while(i<a3):
    inx1=int(redshift3[i]/delta_z)
    inx2=int((M_z3[i]-2)/delta_M)
    dN_dtdzdM3[inx2][inx1]=dN_dtdzdM3[inx2][inx1]+4*np.pi*cMpcyr*W_PS3[i]*X(redshift3[i])**2/delta_M/delta_z
    dN_dt3[inx2][inx1]=dN_dt3[inx2][inx1]+4*np.pi*cMpcyr*W_PS3[i]*X(redshift3[i])**2
    i=i+1;

np.savetxt("dN_dtdzdlogM_PopIII.txt",dN_dtdzdM1,fmt='%f')
np.savetxt("dN_dtdzdlogM_Q3d.txt",dN_dtdzdM2,fmt='%f')
np.savetxt("dN_dtdzdlogM_Q3nod.txt",dN_dtdzdM3,fmt='%f')
'''
redshiftlist=[]
for i in range(0,num):
    for j in range(0,N):
        a=dN_dt2[i][j]=int(dN_dt2[i][j]+0.5)
        #print(a)
        if(a!=0):
            redshiftlist.append(np.random.uniform(j*delta_z,(j+1)*delta_z,a))
np.savetxt("dNdtQ3d.txt",dN_dt2,fmt='%d')
print(redshiftlist)
'''


'''
#generate BBH catalogue
numberPop=(np.sum(dN_dt1))
numberQ3d=(np.sum(dN_dt2))
numberQ3nod=(np.sum(dN_dt3))
print(numberPop,numberQ3d,numberQ3nod)
print(np.sum(dN_dtdzdM1)*delta_z*delta_M,np.sum(dN_dtdzdM2)*delta_z*delta_M,np.sum(dN_dtdzdM3)*delta_z*delta_M)
dN_dt1=dN_dt1/numberPop
dN_dt2=dN_dt2/numberQ3d
dN_dt3=dN_dt3/numberQ3nod
print(np.sum(dN_dt1),np.sum(dN_dt2),np.sum(dN_dt3))
#print(dN_dt1)
#np.savetxt("dN_dt1.txt",dN_dt1,fmt='%.15f')
dNdt1=dN_dt1.flatten()
dNdt2=dN_dt2.flatten()
dNdt3=dN_dt3.flatten()
dNdt1=np.array(dNdt1)
#np.savetxt("dNdt1.txt",dN_dt1,fmt='%.15f')
#print(dNdt1)
print(dNdt1.shape)

index1=np.random.choice(num*N,int(numberPop),p=dNdt1)
index2=np.random.choice(num*N,int(numberQ3d+0.5),p=dNdt2)
index3=np.random.choice(num*N,int(numberQ3nod),p=dNdt3)
#print(index1)
logMredshift1=np.zeros((int(numberPop),2))       #(logM,  redshift)
logMredshift2=np.zeros((int(numberQ3d+0.5),2))
logMredshift3=np.zeros((int(numberQ3nod),2))
i=0
m=10
for j in range(0,m):
    index1=np.random.choice(num*N,int(numberPop),p=dNdt1)
    index2=np.random.choice(num*N,int(numberQ3d+0.5),p=dNdt2)
    index3=np.random.choice(num*N,int(numberQ3nod),p=dNdt3)
    for i in range(len(index1)):
        logMredshift1[i][0]=np.random.uniform(index1[i]//N*delta_M+2,index1[i]//N*delta_M+2+delta_M,1)
        logMredshift1[i][1]=np.random.uniform((index1[i]%N)*delta_z,(index1[i]%N)*delta_z+delta_z,1)
        logMredshift1[i][0]=10**(logMredshift1[i][0])

    for i in range(len(index2)):
        logMredshift2[i][0]=np.random.uniform(index2[i]//N*delta_M+2,index2[i]//N*delta_M+2+delta_M,1)
        logMredshift2[i][1]=np.random.uniform((index2[i]%N)*delta_z,(index2[i]%N)*delta_z+delta_z,1)
        logMredshift2[i][0]=10**(logMredshift2[i][0])

    for i in range(len(index3)):
        logMredshift3[i][0]=np.random.uniform(index3[i]//N*delta_M+2,index3[i]//N*delta_M+2+delta_M,1)
        logMredshift3[i][1]=np.random.uniform((index3[i]%N)*delta_z,(index3[i]%N)*delta_z+delta_z,1)
        logMredshift3[i][0]=10**(logMredshift3[i][0])

logMredshift=np.vstack((logMredshift1,logMredshift2,logMredshift3))
max_inx=np.argmax(logMredshift[:,0])
print('max mass',logMredshift[max_inx][0])
print(logMredshift1.shape,logMredshift2.shape,logMredshift3.shape)
print(logMredshift.shape)
#np.savetxt("Mc10.txt",logMredshift[:,0],fmt='%7.1f')
#np.savetxt("redshift10.txt",logMredshift[:,1],fmt='%2.6f')



#############################
#print(logMredshift1)
plt.scatter(logMredshift1[:,1],logMredshift1[:,0])
plt.ylabel('$\log{[M_c/M_{\odot}]}$')
plt.xlabel('redshift')
plt.show()
plt.scatter(logMredshift2[:,1],logMredshift2[:,0])
plt.ylabel('$\log{[M_c/M_{\odot}]}$')
plt.xlabel('redshift')
plt.show()
plt.scatter(logMredshift3[:,1],logMredshift3[:,0])
plt.ylabel('$\log{[M_c/M_{\odot}]}$')
plt.xlabel('redshift')
plt.show()

plt.hist(logMredshift1[:,1],bins=N,rwidth=0.8,label="PopIII")
plt.xlabel('redshift')
plt.show()
plt.hist(logMredshift2[:,1],bins=N,rwidth=0.8,label="Q3delays")
plt.xlabel('redshift')
plt.show()
plt.hist(logMredshift3[:,1],bins=N,rwidth=0.8,label="Q3nodelays")
plt.xlabel('redshift')
plt.show()
plt.hist(logMredshift1[:,0],bins=num,rwidth=0.8,label="PopIII")
plt.xlabel('$\log{[M_c/M_{\odot}]}$')
plt.show()
plt.hist(logMredshift2[:,0],bins=num,rwidth=0.8,label="Q3delays")
plt.xlabel('$\log{[M_c/M_{\odot}]}$')
plt.show()
plt.hist(logMredshift3[:,0],bins=num,rwidth=0.8,label="Q3nodelays")
plt.xlabel('$\log{[M_c/M_{\odot}]}$')
plt.show()



Mc=np.loadtxt("data/Mc_Q3nod9.txt")
redshift=np.loadtxt("data/redshift_Q3nod9.txt")
for i in range(len(Mc)):
    Mc[i]=math.log(Mc[i],10)
#max_inx=np.argmax(Mc)
#print('max mass',Mc[max_inx])
plt.scatter(redshift,Mc)
plt.ylabel('$[M_c/M_{\odot}]$')
plt.xlabel('redshift')
plt.savefig("data1/scatterQ3nod9.png")
plt.show()
plt.hist(redshift,bins=N,rwidth=0.8)
plt.xlabel('redshift')
plt.savefig("data1/Q3nod9redshift.png")
plt.show()
plt.hist(Mc,bins=num,rwidth=0.8)
plt.xlabel('$[M_c/M_{\odot}]$')
plt.savefig("data1/Q3nod9Mc.png")
plt.show()



np.savetxt("dNdtPopIII.txt",dN_dt1,fmt='%d')
np.savetxt("dNdtQ3d.txt",dN_dt2,fmt='%d')
np.savetxt("dNdtQ3nod.txt",dN_dt3,fmt='%d')
maxinx=np.unravel_index(np.argmax(dN_dtdzdM3),dN_dtdzdM3.shape)
total=dN_dtdzdM3[maxinx[0]][maxinx[1]]**delta_z*delta_M
print(2+(maxinx[0])*delta_M)
while(total<=100/2):
    upper=maxinx[0]+1
    total=total+dN_dtdzdM3[upper][maxinx[1]]*delta_z*delta_M
print(total)
total=dN_dtdzdM3[maxinx[0]][maxinx[1]]*delta_z*delta_M
while(total<=100/2):
    low=maxinx[0]-1
    total=total+dN_dtdzdM3[low][maxinx[1]]*delta_z*delta_M
print(total)
print(low*delta_M,(upper+1)*delta_M)
print(maxinx[1]*delta_z)
#print(maxinx[0]*0.4,maxinx[1],dN_dtdzdM1[maxinx[0]][maxinx[1]])
'''
'''
#heatmap
fig = plt.figure()
axx1 = fig.add_subplot(131)
axx2 = fig.add_subplot(132)
axx3 = fig.add_subplot(133)

#ax.set_xticks(range(len(z)))
plt.gcf().set_size_inches(10, 5)
im1 = axx1.imshow(dN_dtdzdM1,interpolation='bilinear', origin='lower',cmap = plt.cm.hot_r)
im2 = axx2.imshow(dN_dtdzdM2,interpolation='bilinear', origin='lower',cmap = plt.cm.hot_r)
im3 = axx3.imshow(dN_dtdzdM3,interpolation='bilinear', origin='lower',cmap = plt.cm.hot_r)
cbar=fig.colorbar(im1,ax=[axx1,axx2,axx3],orientation="horizontal", pad=0.08)


axx1.set_title('PopIII')
axx2.set_title('Q3delays')
axx3.set_title('Q3nodelays')
#plt.xlim(0,N)
#plt.ylim(0,num)

axx1.set_ylabel('$\log{[M_c/M_{\odot}]}$')
axx1.set_xlabel('redshift')
#axx2.set_ylabel('$\log{[M_c/M_{\odot}]}$')
axx2.set_xlabel('redshift')
#axx3.set_ylabel('$\log{[M_c/M_{\odot}]}$')
axx3.set_xlabel('redshift')

tick=['0','5','10','15','20']
tick1=[0.0/delta_z,5/delta_z,10/delta_z,15/delta_z,20/delta_z]
axx1.set_xticks(tick1)
axx2.set_xticks(tick1)
axx3.set_xticks(tick1)
axx1.set_xticklabels(tick)
axx2.set_xticklabels(tick)
axx3.set_xticklabels(tick)
tick=['2.0','4.25','6.5','8.75','11.0']
tick1=[0,num/4,num/2,num/4*3,num]
axx1.set_yticks(tick1)
axx1.set_yticklabels(tick)
axx2.set_yticks(tick1)
axx2.set_yticklabels(tick)
axx3.set_yticks(tick1)
axx3.set_yticklabels(tick)
plt.show()
'''


#Q3nod
for i in range(0,num):
    for j in range(0,N):
        if(dN_dtdzdM3[i][j]!=0):
            dN_dtdzdM3[i][j]=math.log(dN_dtdzdM3[i][j],10)
        else:
            dN_dtdzdM3[i][j]=-8
'''
#network without lensing
exec('first{}=[227]'.format(0));exec('second{}=[213,210]'.format(0));exec('third{}=[194,277,285,282,223,209,263]'.format(0));exec('four{}=[237,201,225,287,198,232,184,246,269]'.format(0));
exec('first{}=[-1]'.format(1));exec('second{}=[208,194,277]'.format(1));exec('third{}=[290,275,225,238,229,176,216,259,250,283,257]'.format(1));exec('four{}=[213,288,276,192,222,220,272,234]'.format(1));
exec('first{}=[-1]'.format(2));exec('second{}=[221,250]'.format(2));exec('third{}=[209,238,260,282,187,175,193,261,186,178]'.format(2));exec('four{}=[251,196,224,176,195,197,203,258,233,262,290]'.format(2));
exec('first{}=[-1]'.format(3));exec('second{}=[281]'.format(3));exec('third{}=[194,263,249,260,211,245,253,190]'.format(3));exec('four{}=[223,283,227,290,228,191]'.format(3));
exec('first{}=[247,191,196]'.format(4));exec('second{}=[222,226]'.format(4));exec('third{}=[179,212,185,207,210,261,225,218,220,229,254,242,211,224]'.format(4));exec('four{}=[284,186,264,241,215,187,180,275,176,237,199,184,217]'.format(4));
exec('first{}=[-1]'.format(5));exec('second{}=[-1]'.format(5));exec('third{}=[188,198,272,251,219,247,262]'.format(5));exec('four{}=[280,270,223,266,208,281,279,250,286,263,184,257,256,275,177,226]'.format(5));
exec('first{}=[-1]'.format(6));exec('second{}=[218]'.format(6));exec('third{}=[189,197,179,245,285,260,234,267,243,220]'.format(6));exec('four{}=[213,192,238,221,204,181,225,200,266,287]'.format(6));
exec('first{}=[-1]'.format(7));exec('second{}=[256,200]'.format(7));exec('third{}=[278,270,209,288,221,285,215,185,260,189]'.format(7));exec('four{}=[244,232,175,274,226,229,183,230,290,281,234]'.format(7));
exec('first{}=[182]'.format(8));exec('second{}=[246]'.format(8));exec('third{}=[241,237,271,232,197,216,243,176,213,288,263,272]'.format(8));exec('four{}=[178,209,282,228,264,210,183]'.format(8));
exec('first{}=[-1]'.format(9));exec('second{}=[253,247]'.format(9));exec('third{}=[290,225,179,191,189,208,243,271,192,251,274,272,185,262,281,211]'.format(9));exec('four{}=[276,200,255,222,224,242,206,235,254,214,188,190]'.format(9));


datanum=0               
for datanum in range(0,10):
    #read the data                                                     #nine-dimensional matrix(M_c,eta,iota,d_L,alpha,delta,t_c,phi_c,psi)
    redshift1=np.loadtxt("threeyears/redshift%d.txt"%datanum)
    M1=np.loadtxt("threeyears/Mc%d.txt"%datanum)                                 #Msun
    i=0
    while(i<874):
        M1[i]=math.log(M1[i],10)
        i=i+1    
    plt.gcf().set_size_inches(10, 7)
    [X, Y] = np.meshgrid(z, log_M)
    list=np.array([-4,-3,-2.0,-1.0,0.0,0.5,1.0,])
    CQ3nod1=plt.contourf(X, Y, dN_dtdzdM3, list,cmap=plt.cm.Greys,alpha=0.5,extend='min')
    CQ3nod= plt.contour(X, Y, dN_dtdzdM3, list,colors=['k','k','k','k','k','k'])
    plt.clabel(CQ3nod, inline=1, inline_spacing=0,fontsize=8,fmt='%1.1f')

    cbar=plt.colorbar(CQ3nod1)
    cbar.set_label('$\log{[dN/(dtdzdlogM_c)]}$')
    #cbar.set_ticks([-8,-6,-4,-2,0,2])
    #cbar.set_ticklabels(['$0$','$10^{-6}$','$10^{-4}$','$10^{-2}$','$10^0$','$10^2$'])
    a0=plt.scatter(redshift1[525:874],M1[525:874],s=10,color='b')
    exec('first=first{}'.format(datanum))
    exec('second=second{}'.format(datanum))
    exec('third=third{}'.format(datanum))
    exec('four=four{}'.format(datanum))
    print(first,second,third)
    if(third[0]>0):
        a3=plt.scatter(redshift1[[third]],M1[[third]],s=100,color='lightgreen',marker='s',label='1%<σ<5%')     #1%<σ(H0)<5%
    if(second[0]>0):
        a2=plt.scatter(redshift1[[second]],M1[[second]],s=180,color='yellow',edgecolor='orange',marker='*',label='0.5%<σ<1%')                         #0.5%<σ(H0)<1%
    if(first[0]>0):
        a1=plt.scatter(redshift1[[first]],M1[[first]],s=80,color='r',marker='d',label='σ<0.5%')                  #σ(H0)<0.5%
    if(four[0]>0):
        a4=plt.scatter(redshift1[[four]],M1[[four]],s=80,color='w',edgecolor='b',label='$N_host>10^6$')
    
    #plt.grid(True)
    plt.legend(handles=[a1,a2,a3,a4,a0],labels=['$σ_{H_0}<0.5\%$','$0.5\%<σ_{H_0}<1\%$','$1\%<σ_{H_0}<5\%$','$σ_{H_0}>5\%$','$N_{host}>10^6$'])
    plt.title('Q3nod w/o lensing noise in network')
    #plt.title('Q3nod w/o lensing noise in Taiji')
    plt.xlim(0,20)
    #plt.ylim(3.5,10.0)
    plt.ylabel('$\log{[M_c/M_{\odot}]}$', fontsize = 16)
    plt.xlabel('redshift ($z$)', fontsize = 16)
    plt.savefig("threeyears/network/nolens/Q3nod%d.png"%datanum, dpi=300)
    #plt.savefig("data1/taiji/nolens/Q3nod%d.png"%datanum, dpi=300)
    plt.show() 





#taiji without lensing
exec('first{}=[-1]'.format(0));exec('second{}=[-1]'.format(0));exec('third{}=[227]'.format(0));exec('four{}=[-1]'.format(0));
exec('first{}=[-1]'.format(1));exec('second{}=[-1]'.format(1));exec('third{}=[-1]'.format(1));exec('four{}=[-1]'.format(1));
exec('first{}=[-1]'.format(2));exec('second{}=[-1]'.format(2));exec('third{}=[-1]'.format(2));exec('four{}=[260]'.format(2));
exec('first{}=[-1]'.format(3));exec('second{}=[-1]'.format(3));exec('third{}=[-1]'.format(3));exec('four{}=[223]'.format(3));
exec('first{}=[247]'.format(4));exec('second{}=[-1]'.format(4));exec('third{}=[-1]'.format(4));exec('four{}=[-1]'.format(4));
exec('first{}=[-1]'.format(5));exec('second{}=[-1]'.format(5));exec('third{}=[-1]'.format(5));exec('four{}=[-1]'.format(5));
exec('first{}=[-1]'.format(6));exec('second{}=[-1]'.format(6));exec('third{}=[-1]'.format(6));exec('four{}=[-1]'.format(6));
exec('first{}=[-1]'.format(7));exec('second{}=[-1]'.format(7));exec('third{}=[-1]'.format(7));exec('four{}=[-1]'.format(7));
exec('first{}=[182]'.format(8));exec('second{}=[-1]'.format(8));exec('third{}=[-1]'.format(8));exec('four{}=[228]'.format(8));
exec('first{}=[-1]'.format(9));exec('second{}=[-1]'.format(9));exec('third{}=[-1]'.format(9));exec('four{}=[-1]'.format(9));



datanum=0               
for datanum in range(0,10):
    #read the data                                                     #nine-dimensional matrix(M_c,eta,iota,d_L,alpha,delta,t_c,phi_c,psi)
    redshift1=np.loadtxt("data/redshift%d.txt"%datanum)
    M1=np.loadtxt("data/Mc%d.txt"%datanum)                                 #Msun
    i=0
    while(i<291):
        M1[i]=math.log(M1[i],10)
        i=i+1    
    plt.gcf().set_size_inches(10, 7)
    [X, Y] = np.meshgrid(z, log_M)
    list=np.array([-4,-3,-2.0,-1.0,0.0,0.5,1.0,])
    CQ3nod1=plt.contourf(X, Y, dN_dtdzdM3, list,cmap=plt.cm.Greys,alpha=0.5,extend='min')
    CQ3nod= plt.contour(X, Y, dN_dtdzdM3, list,colors=['k','k','k','k','k','k'])
    plt.clabel(CQ3nod, inline=1, inline_spacing=0,fontsize=8,fmt='%1.1f')

    cbar=plt.colorbar(CQ3nod1)
    cbar.set_label('$\log{[dN/(dtdzdlogM_c)]}$')
    #cbar.set_ticks([-8,-6,-4,-2,0,2])
    #cbar.set_ticklabels(['$0$','$10^{-6}$','$10^{-4}$','$10^{-2}$','$10^0$','$10^2$'])
    a0=plt.scatter(redshift1[175:291],M1[175:291],s=10,color='b')
    exec('first=first{}'.format(datanum))
    exec('second=second{}'.format(datanum))
    exec('third=third{}'.format(datanum))
    exec('four=four{}'.format(datanum))
    print(first,second,third)
    if(third[0]>0):
        a3=plt.scatter(redshift1[[third]],M1[[third]],s=100,color='lightgreen',marker='s',label='1%<σ<5%')     #1%<σ(H0)<5%
    if(second[0]>0):
        a2=plt.scatter(redshift1[[second]],M1[[second]],s=180,color='yellow',edgecolor='orange',marker='*',label='0.5%<σ<1%')                         #0.5%<σ(H0)<1%
    if(first[0]>0):
        a1=plt.scatter(redshift1[[first]],M1[[first]],s=80,color='r',marker='d',label='σ<0.5%')                  #σ(H0)<0.5%
    if(four[0]>0):
        a4=plt.scatter(redshift1[[four]],M1[[four]],s=80,color='w',edgecolor='b',label='$N_host>10^6$')
    
    #plt.grid(True)
    plt.legend(handles=[a1,a2,a3,a4,a0],labels=['$σ_{H_0}<0.5\%$','$0.5\%<σ_{H_0}<1\%$','$1\%<σ_{H_0}<5\%$','$σ_{H_0}>5\%$','$N_{host}>10^6$'])
    #plt.title('Q3nod w/o lensing noise in network')
    plt.title('Q3nod w/o lensing noise in Taiji')
    plt.xlim(0,20)
    plt.ylim(3.5,10.0)
    plt.ylabel('$\log{[M_c/M_{\odot}]}$', fontsize = 16)
    plt.xlabel('redshift ($z$)', fontsize = 16)
    #plt.savefig("data/network/nolens/Q3nod%d.png"%datanum, dpi=300)
    plt.savefig("data/taiji/nolens/Q3nod%d.png"%datanum, dpi=300)
    plt.show() 
'''
##########################################################Q3d
for i in range(0,num):
    for j in range(0,N):
        if(dN_dtdzdM2[i][j]!=0):
            dN_dtdzdM2[i][j]=math.log(dN_dtdzdM2[i][j],10)
        else:
            dN_dtdzdM2[i][j]=-8    
'''
#network without lensing 
exec('first{}=[-1]'.format(0));exec('second{}=[-1]'.format(0));exec('third{}=[-1]'.format(0));exec('four{}=[172]'.format(0));
exec('first{}=[-1]'.format(1));exec('second{}=[169,171]'.format(1));exec('third{}=[170,174,173,167]'.format(1));exec('four{}=[-1]'.format(1));
exec('first{}=[-1]'.format(2));exec('second{}=[-1]'.format(2));exec('third{}=[174]'.format(2));exec('four{}=[172]'.format(2));
exec('first{}=[-1]'.format(3));exec('second{}=[-1]'.format(3));exec('third{}=[-1]'.format(3));exec('four{}=[167]'.format(3));
exec('first{}=[-1]'.format(4));exec('second{}=[-1]'.format(4));exec('third{}=[172,168,173]'.format(4));exec('four{}=[-1]'.format(4));
exec('first{}=[-1]'.format(5));exec('second{}=[-1]'.format(5));exec('third{}=[168]'.format(5));exec('four{}=[172]'.format(5));
exec('first{}=[-1]'.format(6));exec('second{}=[-1]'.format(6));exec('third{}=[171]'.format(6));exec('four{}=[170,169]'.format(6));
exec('first{}=[-1]'.format(7));exec('second{}=[-1]'.format(7));exec('third{}=[169]'.format(7));exec('four{}=[170,171]'.format(7));
exec('first{}=[168]'.format(8));exec('second{}=[-1]'.format(8));exec('third{}=[173]'.format(8));exec('four{}=[169]'.format(8));
exec('first{}=[-1]'.format(9));exec('second{}=[171]'.format(9));exec('third{}=[170,167,174]'.format(9));exec('four{}=[172]'.format(9));

#taiji without lensing
exec('first{}=[-1]'.format(0));exec('second{}=[-1]'.format(0));exec('third{}=[-1]'.format(0));exec('four{}=[-1]'.format(0));
exec('first{}=[-1]'.format(1));exec('second{}=[-1]'.format(1));exec('third{}=[-1]'.format(1));exec('four{}=[-1]'.format(1));
exec('first{}=[-1]'.format(2));exec('second{}=[-1]'.format(2));exec('third{}=[-1]'.format(2));exec('four{}=[-1]'.format(2));
exec('first{}=[-1]'.format(3));exec('second{}=[-1]'.format(3));exec('third{}=[-1]'.format(3));exec('four{}=[-1]'.format(3));
exec('first{}=[-1]'.format(4));exec('second{}=[-1]'.format(4));exec('third{}=[-1]'.format(4));exec('four{}=[-1]'.format(4));
exec('first{}=[-1]'.format(5));exec('second{}=[-1]'.format(5));exec('third{}=[-1]'.format(5));exec('four{}=[-1]'.format(5));
exec('first{}=[-1]'.format(6));exec('second{}=[-1]'.format(6));exec('third{}=[-1]'.format(6));exec('four{}=[-1]'.format(6));
exec('first{}=[-1]'.format(7));exec('second{}=[-1]'.format(7));exec('third{}=[-1]'.format(7));exec('four{}=[-1]'.format(7));
exec('first{}=[168]'.format(8));exec('second{}=[-1]'.format(8));exec('third{}=[-1]'.format(8));exec('four{}=[-1]'.format(8));
exec('first{}=[-1]'.format(9));exec('second{}=[-1]'.format(9));exec('third{}=[-1]'.format(9));exec('four{}=[-1]'.format(9));


datanum=0               
for datanum in range(0,10):
    #read the data                                                     #nine-dimensional matrix(M_c,eta,iota,d_L,alpha,delta,t_c,phi_c,psi)
    redshift1=np.loadtxt("data/redshift%d.txt"%datanum)
    M1=np.loadtxt("data/Mc%d.txt"%datanum)                                 #Msun
    i=0
    while(i<291):
        M1[i]=math.log(M1[i],10)
        i=i+1
    plt.gcf().set_size_inches(10, 7)
    [X, Y] = np.meshgrid(z, log_M)
    list=np.array([-4.0,-2.0,-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.1,0.2])
    CQ3d1=plt.contourf(X, Y, dN_dtdzdM2, list,cmap=plt.cm.Greys,alpha=0.5,extend='min')
    CQ3d= plt.contour(X, Y, dN_dtdzdM2, list,colors=['k','k','k','k','k','k'])
    plt.clabel(CQ3d, inline=1, inline_spacing=0,fontsize=8,fmt='%1.1f')
    cbar=plt.colorbar(CQ3d1)
    cbar.set_label('$\log{[dN/(dtdzdlogM_c)]}$')
    #cbar.set_ticks([-8,-6,-4,-2,0,2])
    #cbar.set_ticklabels(['$0$','$10^{-6}$','$10^{-4}$','$10^{-2}$','$10^0$','$10^2$'])
    a0=plt.scatter(redshift1[167:175],M1[167:175],s=10,color='b')
    exec('first=first{}'.format(datanum))
    exec('second=second{}'.format(datanum))
    exec('third=third{}'.format(datanum))
    exec('four=four{}'.format(datanum))
    print(first,second,third)
    if(third[0]>0):
        a3=plt.scatter(redshift1[[third]],M1[[third]],s=100,color='lightgreen',marker='s',label='1%<σ<5%')     #1%<σ(H0)<5%
    if(second[0]>0):
        a2=plt.scatter(redshift1[[second]],M1[[second]],s=180,color='yellow',edgecolor='orange',marker='*',label='0.5%<σ<1%')                         #0.5%<σ(H0)<1%
    if(first[0]>0):
        a1=plt.scatter(redshift1[[first]],M1[[first]],s=80,color='r',marker='d',label='σ<0.5%')                  #σ(H0)<0.5%
    if(four[0]>0):
        a4=plt.scatter(redshift1[[four]],M1[[four]],s=80,color='w',edgecolor='b',label='N_host>10^6')     
    #plt.grid(True)
    plt.legend(handles=[a1,a2,a3,a4,a0],labels=['$σ_{H_0}<0.5\%$','$0.5\%<σ_{H_0}<1\%$','$1\%<σ_{H_0}<5\%$','$σ_{H_0}>5\%$','$N_{host}>10^6$'])
    #plt.title('Q3d w/o lensing noise in network')
    plt.title('Q3d w/o lensing noise in Taiji')
    plt.xlim(0,20)
    plt.ylim(3.0,10.0)
    plt.ylabel('$\log{[M_c/M_{\odot}]}$', fontsize = 16)
    plt.xlabel('redshift ($z$)', fontsize = 16)
    #plt.savefig("data/network/nolens/Q3d%d.png"%datanum,dpi=300)
    plt.savefig("data/taiji/nolens/Q3d%d.png"%datanum,dpi=300)
    plt.show() 
'''
#PopIII
for i in range(0,num):
    for j in range(0,N):
        if(dN_dtdzdM1[i][j]!=0):
            dN_dtdzdM1[i][j]=math.log(dN_dtdzdM1[i][j],10)
        else:
            dN_dtdzdM1[i][j]=-8 
            
'''
#network without lensing
exec('first{}=[-1]'.format(0));exec('second{}=[-1]'.format(0));exec('third{}=[77]'.format(0));exec('four{}=[159]'.format(0));
exec('first{}=[34]'.format(1));exec('second{}=[140]'.format(1));exec('third{}=[-1]'.format(1));exec('four{}=[-1]'.format(1));
exec('first{}=[-1]'.format(2));exec('second{}=[8]'.format(2));exec('third{}=[66,145]'.format(2));exec('four{}=[-1]'.format(2));
exec('first{}=[-1]'.format(3));exec('second{}=[-1]'.format(3));exec('third{}=[-1]'.format(3));exec('four{}=[161,68,63,15]'.format(3));
exec('first{}=[156]'.format(4));exec('second{}=[-1]'.format(4));exec('third{}=[-1]'.format(4));exec('four{}=[54,119]'.format(4));
exec('first{}=[-1]'.format(5));exec('second{}=[-1]'.format(5));exec('third{}=[-1]'.format(5));exec('four{}=[4,0,8]'.format(5));
exec('first{}=[-1]'.format(6));exec('second{}=[-1]'.format(6));exec('third{}=[135]'.format(6));exec('four{}=[8]'.format(6));
exec('first{}=[-1]'.format(7));exec('second{}=[-1]'.format(7));exec('third{}=[-1]'.format(7));exec('four{}=[-1]'.format(7));
exec('first{}=[-1]'.format(8));exec('second{}=[-1]'.format(8));exec('third{}=[9,161,40]'.format(8));exec('four{}=[153]'.format(8));
exec('first{}=[122]'.format(9));exec('second{}=[91]'.format(9));exec('third{}=[-1]'.format(9));exec('four{}=[53]'.format(9));

#Taiji without lensing
exec('first{}=[-1]'.format(0));exec('second{}=[-1]'.format(0));exec('third{}=[-1]'.format(0));exec('four{}=[-1]'.format(0));
exec('first{}=[-1]'.format(1));exec('second{}=[-1]'.format(1));exec('third{}=[-1]'.format(1));exec('four{}=[-1]'.format(1));
exec('first{}=[-1]'.format(2));exec('second{}=[-1]'.format(2));exec('third{}=[-1]'.format(2));exec('four{}=[-1]'.format(2));
exec('first{}=[-1]'.format(3));exec('second{}=[-1]'.format(3));exec('third{}=[-1]'.format(3));exec('four{}=[-1]'.format(3));
exec('first{}=[-1]'.format(4));exec('second{}=[-1]'.format(4));exec('third{}=[156]'.format(4));exec('four{}=[-1]'.format(4));
exec('first{}=[-1]'.format(5));exec('second{}=[-1]'.format(5));exec('third{}=[-1]'.format(5));exec('four{}=[-1]'.format(5));
exec('first{}=[-1]'.format(6));exec('second{}=[-1]'.format(6));exec('third{}=[-1]'.format(6));exec('four{}=[-1]'.format(6));
exec('first{}=[-1]'.format(7));exec('second{}=[-1]'.format(7));exec('third{}=[-1]'.format(7));exec('four{}=[-1]'.format(7));
exec('first{}=[-1]'.format(8));exec('second{}=[-1]'.format(8));exec('third{}=[-1]'.format(8));exec('four{}=[-1]'.format(8));
exec('first{}=[-1]'.format(9));exec('second{}=[-1]'.format(9));exec('third{}=[-1]'.format(9));exec('four{}=[-1]'.format(9));


datanum=0               
for datanum in range(0,10):
    #read the data                                                     #nine-dimensional matrix(M_c,eta,iota,d_L,alpha,delta,t_c,phi_c,psi)
    redshift1=np.loadtxt("data/redshift%d.txt"%datanum)
    M1=np.loadtxt("data/Mc%d.txt"%datanum)                                 #Msun
    i=0
    while(i<291):
        M1[i]=math.log(M1[i],10)
        i=i+1
    plt.scatter(redshift1,M1,s=10,color='r')
    plt.show()
    #等高线
    plt.gcf().set_size_inches(10, 7)
    [X, Y] = np.meshgrid(z, log_M)

    list=np.array([-4,-2.0,-1.0,0.0,1.0,1.5])
    print(list)
    #CPop1=plt.contourf(X, Y, dN_dtdzdM1,list,cmap=plt.cm.Greys,alpha=0.5)
    CPop1=plt.contourf(X, Y, dN_dtdzdM1,list,cmap=plt.cm.Greys,alpha=0.5,extend='min')
    CPop= plt.contour(X, Y, dN_dtdzdM1, list,colors=['k','k','k','k','k','k'])
    plt.clabel(CPop, inline=1, inline_spacing=0,fontsize=8,fmt='%1.1f')
    cbar=plt.colorbar(CPop1)
    cbar.set_label('$\log{[dN/(dtdzdlogM_c)]}$')
    #cbar.set_ticks([-8,-6,-4,-2,0,2])
    #cbar.set_ticklabels(['0','$10^{-4}$','$10^{-2}$','$10^0$','$10^2$'])
    a0=plt.scatter(redshift1[0:166],M1[0:166],s=10,color='b')
    exec('first=first{}'.format(datanum))
    exec('second=second{}'.format(datanum))
    exec('third=third{}'.format(datanum))
    exec('four=four{}'.format(datanum))
    print(first,second,third)
    if(third[0]>0):
        a3=plt.scatter(redshift1[[third]],M1[[third]],s=100,color='lightgreen',marker='s',label='1%<σ<5%')     #1%<σ(H0)<5%
    if(second[0]>0):
        a2=plt.scatter(redshift1[[second]],M1[[second]],s=180,color='yellow',edgecolor='orange',marker='*',label='0.5%<σ<1%')                         #0.5%<σ(H0)<1%
    if(first[0]>0):
        a1=plt.scatter(redshift1[[first]],M1[[first]],s=80,color='r',marker='d',label='σ<0.5%')                  #σ(H0)<0.5%
    if(four[0]>0):
        a4=plt.scatter(redshift1[[four]],M1[[four]],s=80,color='w',edgecolor='b',label='N_host>10^6')      
    #plt.grid(True)
    plt.legend(handles=[a1,a2,a3,a4,a0],labels=['$σ_{H_0}<0.5\%$','$0.5\%<σ_{H_0}<1\%$','$1\%<σ_{H_0}<5\%$','$σ_{H_0}>5\%$','$N_{host}>10^6$'])
    #plt.title('PopIII w/o lensing noise in network')
    plt.title('PopIII w/o lensing noise in Taiji')
    plt.xlim(0,20)
    plt.ylim(2.3,10.0)
    plt.ylabel('$\log{[M_c/M_{\odot}]}$', fontsize = 16)
    plt.xlabel('redshift ($z$)', fontsize = 16)
    #plt.savefig("data/network/nolens/PopIII%d.png"%datanum, dpi=300)
    plt.savefig("data/taiji/nolens/PopIII%d.png"%datanum, dpi=300)
    plt.show() 
'''





#lensing########################

###########Q3nod
'''
for i in range(0,num):
    for j in range(0,N):
        if(dN_dtdzdM3[i][j]!=0):
            dN_dtdzdM3[i][j]=math.log(dN_dtdzdM3[i][j],10)
        else:
            dN_dtdzdM3[i][j]=-8
'''



#network with lensing  in 3 years
exec('first{}=[-1]'.format(0));exec('second{}=[-1]'.format(0));exec('third{}=[709,737,601]'.format(0));exec('four{}=[864,686,584,770,861,660]'.format(0));
exec('first{}=[-1]'.format(1));exec('second{}=[-1]'.format(1));exec('third{}=[-1]'.format(1));exec('four{}=[764,819,687,555,808,790,834]'.format(1));
exec('first{}=[617,540]'.format(2));exec('second{}=[-1]'.format(2));exec('third{}=[614]'.format(2));exec('four{}=[546,787,627,800,696]'.format(2));
exec('first{}=[-1]'.format(3));exec('second{}=[-1]'.format(3));exec('third{}=[578,675,867]'.format(3));exec('four{}=[602,620,733,638,741]'.format(3));
exec('first{}=[-1]'.format(4));exec('second{}=[-1]'.format(4));exec('third{}=[599,631]'.format(4));exec('four{}=[750,659,799,582,835]'.format(4));
exec('first{}=[-1]'.format(5));exec('second{}=[-1]'.format(5));exec('third{}=[-1]'.format(5));exec('four{}=[664,676,847]'.format(5));
exec('first{}=[-1]'.format(6));exec('second{}=[-1]'.format(6));exec('third{}=[648,788]'.format(6));exec('four{}=[558,695,832,612,864]'.format(6));
exec('first{}=[-1]'.format(7));exec('second{}=[-1]'.format(7));exec('third{}=[596,581,706,640,849,748]'.format(7));exec('four{}=[691,754,817,650]'.format(7));
exec('first{}=[-1]'.format(8));exec('second{}=[-1]'.format(8));exec('third{}=[721,767]'.format(8));exec('four{}=[580,745,679,665,823,552,609,638,737]'.format(8));
exec('first{}=[-1]'.format(9));exec('second{}=[-1]'.format(9));exec('third{}=[552,85]'.format(9));exec('four{}=[621,549,771]'.format(9));
'''
#taiji with lensing in 3 years
exec('first{}=[-1]'.format(0));exec('second{}=[-1]'.format(0));exec('third{}=[709]'.format(0));exec('four{}=[861,864]'.format(0));
exec('first{}=[-1]'.format(1));exec('second{}=[-1]'.format(1));exec('third{}=[-1]'.format(1));exec('four{}=[-1]'.format(1));
exec('first{}=[617,540]'.format(2));exec('second{}=[-1]'.format(2));exec('third{}=[614]'.format(2));exec('four{}=[-1]'.format(2));
exec('first{}=[-1]'.format(3));exec('second{}=[-1]'.format(3));exec('third{}=[-1]'.format(3));exec('four{}=[-1]'.format(3));
exec('first{}=[-1]'.format(4));exec('second{}=[-1]'.format(4));exec('third{}=[-1]'.format(4));exec('four{}=[750]'.format(4));
exec('first{}=[-1]'.format(5));exec('second{}=[-1]'.format(5));exec('third{}=[-1]'.format(5));exec('four{}=[-1]'.format(5));
exec('first{}=[-1]'.format(6));exec('second{}=[-1]'.format(6));exec('third{}=[-1]'.format(6));exec('four{}=[558,648]'.format(6));
exec('first{}=[-1]'.format(7));exec('second{}=[-1]'.format(7));exec('third{}=[596]'.format(7));exec('four{}=[748,849]'.format(7));
exec('first{}=[-1]'.format(8));exec('second{}=[-1]'.format(8));exec('third{}=[-1]'.format(8));exec('four{}=[552]'.format(8));
exec('first{}=[-1]'.format(9));exec('second{}=[-1]'.format(9));exec('third{}=[-1]'.format(9));exec('four{}=[-1]'.format(9));
'''
datanum=0               
for datanum in range(0,10):
    #read the data                                                     #nine-dimensional matrix(M_c,eta,iota,d_L,alpha,delta,t_c,phi_c,psi)
    redshift1=np.loadtxt("threeyears/redshift%d.txt"%datanum)
    M1=np.loadtxt("threeyears/Mc%d.txt"%datanum)                                 #Msun
    i=0
    while(i<874):
        M1[i]=math.log(M1[i],10)
        i=i+1    
    plt.gcf().set_size_inches(10, 7)
    [X, Y] = np.meshgrid(z, log_M)
    list=np.array([-4,-3,-2.0,-1.0,0.0,0.5,1.0,])
    CQ3nod1=plt.contourf(X, Y, dN_dtdzdM3, list,cmap=plt.cm.Greys,alpha=0.5,extend='min')
    CQ3nod= plt.contour(X, Y, dN_dtdzdM3, list,colors=['k','k','k','k','k','k'])
    plt.clabel(CQ3nod, inline=1, inline_spacing=0,fontsize=8,fmt='%1.1f')

    cbar=plt.colorbar(CQ3nod1)
    cbar.set_label('$\log{[dN/(dtdzdlogM_c)]}$')
    #cbar.set_ticks([-8,-6,-4,-2,0,2])
    #cbar.set_ticklabels(['$0$','$10^{-6}$','$10^{-4}$','$10^{-2}$','$10^0$','$10^2$'])
    a0=plt.scatter(redshift1[525:874],M1[525:874],s=10,color='b')
    exec('first=first{}'.format(datanum))
    exec('second=second{}'.format(datanum))
    exec('third=third{}'.format(datanum))
    exec('four=four{}'.format(datanum))
    print(first,second,third)
    if(third[0]>0):
        a3=plt.scatter(redshift1[[third]],M1[[third]],s=100,color='lightgreen',marker='s',label='1%<σ<5%')     #1%<σ(H0)<5%
    if(second[0]>0):
        a2=plt.scatter(redshift1[[second]],M1[[second]],s=180,color='yellow',edgecolor='orange',marker='*',label='0.5%<σ<1%')                         #0.5%<σ(H0)<1%
    if(first[0]>0):
        a1=plt.scatter(redshift1[[first]],M1[[first]],s=80,color='r',marker='d',label='σ<0.5%')                  #σ(H0)<0.5%    
    if(four[0]>0):
        a4=plt.scatter(redshift1[[four]],M1[[four]],s=80,color='w',edgecolor='b',label='σ>5%')
    #plt.grid(True)
    plt.legend(handles=[a1,a2,a3,a4,a0],labels=['$σ_{H_0}<0.5\%$','$0.5\%<σ_{H_0}<1\%$','$1\%<σ_{H_0}<5\%$','$σ_{H_0}>5\%$','$N_{host}>10^6$'])
    #plt.legend(loc = 'upper right')
    plt.title('Q3nod w/ lensing noise in network')
    #plt.title('Q3nod w/ lensing noise in Taiji')
    plt.xlim(0,20)
    plt.ylim(3.5,10.0)
    plt.ylabel('$\log{[M_c/M_{\odot}]}$', fontsize = 16)
    plt.xlabel('redshift ($z$)', fontsize = 16)
    plt.savefig("threeyears/network/lens/Q3nodlens%d.png"%datanum, dpi=300)
    #plt.savefig("threeyears/taiji/lens/Q3nodlens%d.png"%datanum, dpi=300)
    plt.show() 


#############PopIII
'''    
for i in range(0,num):
    for j in range(0,N):
        if(dN_dtdzdM1[i][j]!=0):
            dN_dtdzdM1[i][j]=math.log(dN_dtdzdM1[i][j],10)
        else:
            dN_dtdzdM1[i][j]=-8
'''

#network with lensing
exec('first{}=[-1]'.format(0));exec('second{}=[-1]'.format(0));exec('third{}=[149]'.format(0));exec('four{}=[-1]'.format(0));
exec('first{}=[-1]'.format(1));exec('second{}=[-1]'.format(1));exec('third{}=[360]'.format(1));exec('four{}=[-1]'.format(1));
exec('first{}=[-1]'.format(2));exec('second{}=[-1]'.format(2));exec('third{}=[-1]'.format(2));exec('four{}=[316,22,80,315]'.format(2));
exec('first{}=[-1]'.format(3));exec('second{}=[-1]'.format(3));exec('third{}=[-1]'.format(3));exec('four{}=[-1]'.format(3));
exec('first{}=[238]'.format(4));exec('second{}=[-1]'.format(4));exec('third{}=[171]'.format(4));exec('four{}=[211,29,5]'.format(4));
exec('first{}=[-1]'.format(5));exec('second{}=[-1]'.format(5));exec('third{}=[167,292]'.format(5));exec('four{}=[230,59]'.format(5));
exec('first{}=[-1]'.format(6));exec('second{}=[-1]'.format(6));exec('third{}=[42]'.format(6));exec('four{}=[411]'.format(6));
exec('first{}=[-1]'.format(7));exec('second{}=[-1]'.format(7));exec('third{}=[-1]'.format(7));exec('four{}=[126,291]'.format(7));
exec('first{}=[-1]'.format(8));exec('second{}=[-1]'.format(8));exec('third{}=[-1]'.format(8));exec('four{}=[494]'.format(8));
exec('first{}=[23]'.format(9));exec('second{}=[-1]'.format(9));exec('third{}=[-1]'.format(9));exec('four{}=[-1]'.format(9));
'''

#Taiji with lensing
exec('first{}=[-1]'.format(0));exec('second{}=[-1]'.format(0));exec('third{}=[-1]'.format(0));exec('four{}=[-1]'.format(0));
exec('first{}=[-1]'.format(1));exec('second{}=[-1]'.format(1));exec('third{}=[-1]'.format(1));exec('four{}=[-1]'.format(1));
exec('first{}=[-1]'.format(2));exec('second{}=[-1]'.format(2));exec('third{}=[-1]'.format(2));exec('four{}=[-1]'.format(2));
exec('first{}=[-1]'.format(3));exec('second{}=[-1]'.format(3));exec('third{}=[-1]'.format(3));exec('four{}=[-1]'.format(3));
exec('first{}=[-1]'.format(4));exec('second{}=[238]'.format(4));exec('third{}=[-1]'.format(4));exec('four{}=[29]'.format(4));
exec('first{}=[-1]'.format(5));exec('second{}=[-1]'.format(5));exec('third{}=[167]'.format(5));exec('four{}=[-1]'.format(5));
exec('first{}=[-1]'.format(6));exec('second{}=[-1]'.format(6));exec('third{}=[-1]'.format(6));exec('four{}=[42]'.format(6));
exec('first{}=[-1]'.format(7));exec('second{}=[-1]'.format(7));exec('third{}=[-1]'.format(7));exec('four{}=[-1]'.format(7));
exec('first{}=[-1]'.format(8));exec('second{}=[-1]'.format(8));exec('third{}=[-1]'.format(8));exec('four{}=[-1]'.format(8));
exec('first{}=[-1]'.format(9));exec('second{}=[23]'.format(9));exec('third{}=[-1]'.format(9));exec('four{}=[-1]'.format(9));
'''
datanum=0               
for datanum in range(0,0):
    #read the data                                                     #nine-dimensional matrix(M_c,eta,iota,d_L,alpha,delta,t_c,phi_c,psi)
    redshift1=np.loadtxt("threeyears/redshift%d.txt"%datanum)
    M1=np.loadtxt("threeyears/Mc%d.txt"%datanum)                                 #Msun
    i=0
    while(i<874):
        M1[i]=math.log(M1[i],10)
        i=i+1
    plt.scatter(redshift1,M1,s=10,color='r')
    plt.show()
    #等高线
    plt.gcf().set_size_inches(10, 7)
    [X, Y] = np.meshgrid(z, log_M)

    list=np.array([-4,-2.0,-1.0,0.0,1.0,1.5])
    #print(list)
    #CPop1=plt.contourf(X, Y, dN_dtdzdM1,list,cmap=plt.cm.Greys,alpha=0.5)
    CPop1=plt.contourf(X, Y, dN_dtdzdM1,list,cmap=plt.cm.Greys,alpha=0.5,extend='min')
    CPop= plt.contour(X, Y, dN_dtdzdM1, list,colors=['k','k','k','k','k','k'])
    plt.clabel(CPop, inline=1, inline_spacing=0,fontsize=8,fmt='%1.1f')
    cbar=plt.colorbar(CPop1)
    cbar.set_label('$\log{[dN/(dtdzdlogM_c)]}$')
    #cbar.set_ticks([-8,-6,-4,-2,0,2])
    #cbar.set_ticklabels(['0','$10^{-4}$','$10^{-2}$','$10^0$','$10^2$'])
    a0=plt.scatter(redshift1[0:501],M1[0:501],s=10,color='b')
    exec('first=first{}'.format(datanum))
    exec('second=second{}'.format(datanum))
    exec('third=third{}'.format(datanum))
    exec('four=four{}'.format(datanum))
    print(first,second,third)
    if(third[0]>0):
        a3=plt.scatter(redshift1[[third]],M1[[third]],s=100,color='lightgreen',marker='s',label='1%<σ<5%')     #1%<σ(H0)<5%
    if(second[0]>0):
        a2=plt.scatter(redshift1[[second]],M1[[second]],s=180,color='yellow',edgecolor='orange',marker='*',label='0.5%<σ<1%')                         #0.5%<σ(H0)<1%
    if(first[0]>0):
        a1=plt.scatter(redshift1[[first]],M1[[first]],s=80,color='r',marker='d',label='σ<0.5%')                  #σ(H0)<0.5%
    if(four[0]>0):
        a4=plt.scatter(redshift1[[four]],M1[[four]],s=80,color='w',edgecolor='b',label='σ>5%')
    #plt.grid(True)
    plt.legend(handles=[a1,a2,a3,a4,a0],labels=['$σ_{H_0}<0.5\%$','$0.5\%<σ_{H_0}<1\%$','$1\%<σ_{H_0}<5\%$','$σ_{H_0}>5\%$','$N_{host}>10^6$'])
    plt.title('PopIII w/ lensing noise in network')
    #plt.title('PopIII w/ lensing noise in Taiji')
    plt.xlim(0,20)
    plt.ylim(2.3,10.0)
    plt.ylabel('$\log{[M_c/M_{\odot}]}$', fontsize = 16)
    plt.xlabel('redshift ($z$)', fontsize = 16)
    plt.savefig("threeyears/network/lens/PopIIIlens%d.png"%datanum, dpi=300)
    #plt.savefig("threeyears/taiji/lens/PopIIIlens%d.png"%datanum, dpi=300)
    plt.show() 



############Q3d
'''
for i in range(0,num):
    for j in range(0,N):
        if(dN_dtdzdM2[i][j]!=0):
            dN_dtdzdM2[i][j]=math.log(dN_dtdzdM2[i][j],10)
        else:
            dN_dtdzdM2[i][j]=-8    
'''

#network wirh lensing
exec('first{}=[-1]'.format(0));exec('second{}=[-1]'.format(0));exec('third{}=[-1]'.format(0));exec('four{}=[509]'.format(0));
exec('first{}=[-1]'.format(1));exec('second{}=[-1]'.format(1));exec('third{}=[510]'.format(1));exec('four{}=[-1]'.format(1));
exec('first{}=[-1]'.format(2));exec('second{}=[-1]'.format(2));exec('third{}=[-1]'.format(2));exec('four{}=[503]'.format(2));
exec('first{}=[-1]'.format(3));exec('second{}=[-1]'.format(3));exec('third{}=[-1]'.format(3));exec('four{}=[-1]'.format(3));
exec('first{}=[-1]'.format(4));exec('second{}=[-1]'.format(4));exec('third{}=[511,522]'.format(4));exec('four{}=[-1]'.format(4));
exec('first{}=[-1]'.format(5));exec('second{}=[-1]'.format(5));exec('third{}=[-1]'.format(5));exec('four{}=[502,513]'.format(5));
exec('first{}=[-1]'.format(6));exec('second{}=[-1]'.format(6));exec('third{}=[-1]'.format(6));exec('four{}=[-1]'.format(6));
exec('first{}=[-1]'.format(7));exec('second{}=[-1]'.format(7));exec('third{}=[-1]'.format(7));exec('four{}=[521]'.format(7));
exec('first{}=[-1]'.format(8));exec('second{}=[-1]'.format(8));exec('third{}=[523]'.format(8));exec('four{}=[503,519]'.format(8));
exec('first{}=[-1]'.format(9));exec('second{}=[-1]'.format(9));exec('third{}=[-1]'.format(9));exec('four{}=[519]'.format(9));
'''
#taiji with lensing
exec('first{}=[-1]'.format(0));exec('second{}=[-1]'.format(0));exec('third{}=[-1]'.format(0));exec('four{}=[-1]'.format(0));
exec('first{}=[-1]'.format(1));exec('second{}=[-1]'.format(1));exec('third{}=[-1]'.format(1));exec('four{}=[-1]'.format(1));
exec('first{}=[-1]'.format(2));exec('second{}=[-1]'.format(2));exec('third{}=[-1]'.format(2));exec('four{}=[503]'.format(2));
exec('first{}=[-1]'.format(3));exec('second{}=[-1]'.format(3));exec('third{}=[-1]'.format(3));exec('four{}=[-1]'.format(3));
exec('first{}=[-1]'.format(4));exec('second{}=[-1]'.format(4));exec('third{}=[511]'.format(4));exec('four{}=[-1]'.format(4));
exec('first{}=[-1]'.format(5));exec('second{}=[-1]'.format(5));exec('third{}=[-1]'.format(5));exec('four{}=[-1]'.format(5));
exec('first{}=[-1]'.format(6));exec('second{}=[-1]'.format(6));exec('third{}=[-1]'.format(6));exec('four{}=[-1]'.format(6));
exec('first{}=[-1]'.format(7));exec('second{}=[-1]'.format(7));exec('third{}=[-1]'.format(7));exec('four{}=[-1]'.format(7));
exec('first{}=[-1]'.format(8));exec('second{}=[-1]'.format(8));exec('third{}=[-1]'.format(8));exec('four{}=[-1]'.format(8));
exec('first{}=[-1]'.format(9));exec('second{}=[-1]'.format(9));exec('third{}=[-1]'.format(9));exec('four{}=[-1]'.format(9));
'''
datanum=0               
for datanum in range(0,0):
    #read the data                                                     #nine-dimensional matrix(M_c,eta,iota,d_L,alpha,delta,t_c,phi_c,psi)
    redshift1=np.loadtxt("threeyears/redshift%d.txt"%datanum)
    M1=np.loadtxt("threeyears/Mc%d.txt"%datanum)                                 #Msun
    i=0
    while(i<874):
        M1[i]=math.log(M1[i],10)
        i=i+1
    plt.gcf().set_size_inches(10, 7)
    [X, Y] = np.meshgrid(z, log_M)
    list=np.array([-4.0,-2.0,-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.1,0.2])
    CQ3d1=plt.contourf(X, Y, dN_dtdzdM2, list,cmap=plt.cm.Greys,alpha=0.5,extend='min')
    CQ3d= plt.contour(X, Y, dN_dtdzdM2, list,colors=['k','k','k','k','k','k'])
    plt.clabel(CQ3d, inline=1, inline_spacing=0,fontsize=8,fmt='%1.1f')
    cbar=plt.colorbar(CQ3d1)
    cbar.set_label('$\log{[dN/(dtdzdlogM_c)]}$')
    #cbar.set_ticks([-8,-6,-4,-2,0,2])
    #cbar.set_ticklabels(['$0$','$10^{-6}$','$10^{-4}$','$10^{-2}$','$10^0$','$10^2$'])
    a0=plt.scatter(redshift1[501:525],M1[501:525],s=10,color='b')
    exec('first=first{}'.format(datanum))
    exec('second=second{}'.format(datanum))
    exec('third=third{}'.format(datanum))
    exec('four=four{}'.format(datanum))
    print(first,second,third)
    if(third[0]>0):
        a3=plt.scatter(redshift1[[third]],M1[[third]],s=100,color='lightgreen',marker='s',label='1%<σ<5%')     #1%<σ(H0)<5%
    if(second[0]>0):
        a2=plt.scatter(redshift1[[second]],M1[[second]],s=180,color='yellow',edgecolor='orange',marker='*',label='0.5%<σ<1%')                         #0.5%<σ(H0)<1%
    if(first[0]>0):
        a1=plt.scatter(redshift1[[first]],M1[[first]],s=80,color='r',marker='d',label='σ<0.5%')                  #σ(H0)<0.5%
    if(four[0]>0):
        a4=plt.scatter(redshift1[[four]],M1[[four]],s=80,color='w',edgecolor='b',label='σ>5%')
    #plt.grid(True)
    plt.legend(handles=[a1,a2,a3,a4,a0],labels=['$σ_{H_0}<0.5\%$','$0.5\%<σ_{H_0}<1\%$','$1\%<σ_{H_0}<5\%$','$σ_{H_0}>5\%$','$N_{host}>10^6$'])
    plt.title('Q3d w/ lensing noise in network')
    #plt.title('Q3d w/ lensing noise in Taiji')
    plt.xlim(0,20)
    plt.ylim(3.0,10.0)
    plt.ylabel('$\log{[M_c/M_{\odot}]}$', fontsize = 16)
    plt.xlabel('redshift ($z$)', fontsize = 16)
    plt.savefig("threeyears/network/lens/Q3dlens%d.png"%datanum,dpi=300 )
    #plt.savefig("threeyears/taiji/lens/Q3dlens%d.png"%datanum,dpi=300 )
    plt.show() 



##########lensing   in  five years
############Q3nod
'''
for i in range(0,num):
    for j in range(0,N):
        if(dN_dtdzdM3[i][j]!=0):
            dN_dtdzdM3[i][j]=math.log(dN_dtdzdM3[i][j],10)
        else:
            dN_dtdzdM3[i][j]=-8
'''




#network with lensing in  5years
exec('first{}=[1307,1106]'.format(0));exec('second{}=[0]'.format(0));exec('third{}=[1372,1415,1032,1147]'.format(0));exec('four{}=[1159,951,978,1414,1013,933]'.format(0));
exec('first{}=[-1]'.format(1));exec('second{}=[-1]'.format(1));exec('third{}=[1083]'.format(1));exec('four{}=[1248,1366,967,899,895,983]'.format(1));
exec('first{}=[1343]'.format(2));exec('second{}=[-1]'.format(2));exec('third{}=[921,1064,1047]'.format(2));exec('four{}=[1216,889,1442,958,988,1110]'.format(2));
exec('first{}=[-1]'.format(3));exec('second{}=[1184]'.format(3));exec('third{}=[1097,1106,1164,1032,1450]'.format(3));exec('four{}=[1061,1256,1232,1370,911,1272,1036,1276,1107,1373]'.format(3));
exec('first{}=[1181]'.format(4));exec('second{}=[1016]'.format(4));exec('third{}=[1287,1197,1408]'.format(4));exec('four{}=[992,971,1445,1203]'.format(4));
exec('first{}=[902]'.format(5));exec('second{}=[1024]'.format(5));exec('third{}=[1298,1293,1031,1442,992]'.format(5));exec('four{}=[929,883,1223,939,1050,1102,959]'.format(5));
exec('first{}=[-1]'.format(6));exec('second{}=[-1]'.format(6));exec('third{}=[1428,898,1453,1219]'.format(6));exec('four{}=[993,1370,1355,1393,1400,1345]'.format(6));
exec('first{}=[-1]'.format(7));exec('second{}=[-1]'.format(7));exec('third{}=[1271,1309,1257,1338,1006,1417,1196]'.format(7));exec('four{}=[987,883,993,1028,894,1292,1287,884,1439,1047,1245]'.format(7));
exec('first{}=[-1]'.format(8));exec('second{}=[-1]'.format(8));exec('third{}=[1037,1272]'.format(8));exec('four{}=[1044,1082,875,888,912,941,1086,1265,982,1053,1288,979]'.format(8));
exec('first{}=[-1]'.format(9));exec('second{}=[-1]'.format(9));exec('third{}=[1341,1141]'.format(9));exec('four{}=[1074,1354,1029]'.format(9));

datanum=0               
for datanum in range(0,5):
    #read the data                                                     #nine-dimensional matrix(M_c,eta,iota,d_L,alpha,delta,t_c,phi_c,psi)
    redshift1=np.loadtxt("fiveyears/redshift%d.txt"%datanum)
    M1=np.loadtxt("fiveyears/Mc%d.txt"%datanum)                                 #Msun
    i=0
    while(i<1456):
        M1[i]=math.log(M1[i],10)
        i=i+1    
    sns.set_style("white")
    #plt.figure(figsize=(12,10))
    plt.gcf().set_size_inches(10, 7)
    [X, Y] = np.meshgrid(z, log_M)
    list=np.array([-4,-3,-2.0,-1.0,0.0,0.5,1.0,])
    CQ3nod1=plt.contourf(X, Y, dN_dtdzdM3, list,cmap=plt.cm.Greys,alpha=0.5,extend='min')
    CQ3nod= plt.contour(X, Y, dN_dtdzdM3, list,colors=['k','k','k','k','k','k'])
    plt.clabel(CQ3nod, inline=1, inline_spacing=0,fontsize=10,fmt='%1.1f')

    cbar=plt.colorbar(CQ3nod1)
    cbar.set_label('$\log{[dN/(dtdzdlogM_c)]}$',fontsize=20)
    cbar.ax.tick_params(labelsize=15)
    #cbar.set_ticks([-8,-6,-4,-2,0,2])
    #cbar.set_ticklabels(['$0$','$10^{-6}$','$10^{-4}$','$10^{-2}$','$10^0$','$10^2$'])
    a0=plt.scatter(redshift1[875:1456],M1[875:1456],s=10,color='b')
    exec('first=first{}'.format(datanum))
    exec('second=second{}'.format(datanum))
    exec('third=third{}'.format(datanum))
    exec('four=four{}'.format(datanum))
    print(first,second,third)
    if(third[0]>0):
        a3=plt.scatter(redshift1[[third]],M1[[third]],s=100,color='lightgreen',marker='s',label='1%<σ<5%')     #1%<σ(H0)<5%
    if(second[0]>0):
        a2=plt.scatter(redshift1[[second]],M1[[second]],s=180,color='yellow',edgecolor='orange',marker='*',label='0.5%<σ<1%')                         #0.5%<σ(H0)<1%
    if(first[0]>0):
        a1=plt.scatter(redshift1[[first]],M1[[first]],s=80,color='r',marker='d',label='σ<0.5%')                  #σ(H0)<0.5%    
    if(four[0]>0):
        a4=plt.scatter(redshift1[[four]],M1[[four]],s=80,color='w',edgecolor='b',label='σ>5%')
    #plt.grid(True)
    font={'size':18}    
    plt.legend(handles=[a1,a2,a3,a4,a0],labels=['$σ_{H_0}<0.5\%$','$0.5\%<σ_{H_0}<1\%$','$1\%<σ_{H_0}<5\%$','$σ_{H_0}>5\%$','$N_{host}>10^6$'],prop=font)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
   
    plt.title('Q3nod in network',fontdict={'size': 30})
    #plt.title('Q3nod w/ lensing noise in Taiji')
    plt.xlim(0,20)
    plt.ylim(3.5,10.0)
    plt.ylabel('$\log_{10}{[M_c/M_{\odot}]}$', fontsize = 30)
    plt.xlabel('Redshift ($z$)', fontsize =30)
    plt.savefig("fiveyears/network/lens/network_Q3nod_lens%d.pdf"%datanum)
    #plt.savefig("fiveyears/taiji/lens/Q3nodlens%d.png"%datanum, dpi=300)
    plt.show() 


###########PopIII
'''    
for i in range(0,num):
    for j in range(0,N):
        if(dN_dtdzdM1[i][j]!=0):
            dN_dtdzdM1[i][j]=math.log(dN_dtdzdM1[i][j],10)
        else:
            dN_dtdzdM1[i][j]=-8
'''

#network with lensing in five years
exec('first{}=[-1]'.format(0));exec('second{}=[-1]'.format(0));exec('third{}=[58]'.format(0));exec('four{}=[758,145,632,662,305]'.format(0));
exec('first{}=[-1]'.format(1));exec('second{}=[-1]'.format(1));exec('third{}=[-1]'.format(1));exec('four{}=[308,388,774,791]'.format(1));
exec('first{}=[316]'.format(2));exec('second{}=[-1]'.format(2));exec('third{}=[800]'.format(2));exec('four{}=[477,428,737,33]'.format(2));
exec('first{}=[-1]'.format(3));exec('second{}=[-1]'.format(3));exec('third{}=[197]'.format(3));exec('four{}=[217]'.format(3));
exec('first{}=[590]'.format(4));exec('second{}=[-1]'.format(4));exec('third{}=[222]'.format(4));exec('four{}=[142,446,532]'.format(4));
exec('first{}=[-1]'.format(5));exec('second{}=[-1]'.format(5));exec('third{}=[-1]'.format(5));exec('four{}=[358]'.format(5));
exec('first{}=[-1]'.format(6));exec('second{}=[-1]'.format(6));exec('third{}=[555]'.format(6));exec('four{}=[27]'.format(6));
exec('first{}=[517]'.format(7));exec('second{}=[-1]'.format(7));exec('third{}=[-1]'.format(7));exec('four{}=[718,119]'.format(7));
exec('first{}=[574]'.format(8));exec('second{}=[-1]'.format(8));exec('third{}=[-1]'.format(8));exec('four{}=[279,565]'.format(8));
exec('first{}=[-1]'.format(9));exec('second{}=[-1]'.format(9));exec('third{}=[194]'.format(9));exec('four{}=[329,488]'.format(9));

datanum=0               
for datanum in range(4,5):
    #read the data                                                     #nine-dimensional matrix(M_c,eta,iota,d_L,alpha,delta,t_c,phi_c,psi)
    redshift1=np.loadtxt("fiveyears/redshift%d.txt"%datanum)
    M1=np.loadtxt("fiveyears/Mc%d.txt"%datanum)                                 #Msun
    i=0
    while(i<1456):
        M1[i]=math.log(M1[i],10)
        i=i+1
    #plt.scatter(redshift1,M1,s=10,color='r')
    #plt.show()
    #等高线
    sns.set_style("white")
    #plt.figure(figsize=(12,10))
    plt.gcf().set_size_inches(10, 7)
    [X, Y] = np.meshgrid(z, log_M)

    list=np.array([-4,-2.0,-1.0,0.0,1.0,1.5])
    #print(list)
    #CPop1=plt.contourf(X, Y, dN_dtdzdM1,list,cmap=plt.cm.Greys,alpha=0.5)
    CPop1=plt.contourf(X, Y, dN_dtdzdM1,list,cmap=plt.cm.Greys,alpha=0.5,extend='min')
    CPop= plt.contour(X, Y, dN_dtdzdM1, list,colors=['k','k','k','k','k','k'])
    plt.clabel(CPop, inline=1, inline_spacing=0,fontsize=11,fmt='%1.1f')
    cbar=plt.colorbar(CPop1)
    cbar.set_label('$\log{[dN/(dtdzdlogM_c)]}$',fontsize=20)
    cbar.ax.tick_params(labelsize=15)
    #cbar.set_ticks([-8,-6,-4,-2,0,2])
    #cbar.set_ticklabels(['0','$10^{-4}$','$10^{-2}$','$10^0$','$10^2$'])
    a0=plt.scatter(redshift1[0:835],M1[0:835],s=10,color='b')
    exec('first=first{}'.format(datanum))
    exec('second=second{}'.format(datanum))
    exec('third=third{}'.format(datanum))
    exec('four=four{}'.format(datanum))
    print(first,second,third)
    if(third[0]>0):
        a3=plt.scatter(redshift1[[third]],M1[[third]],s=100,color='lightgreen',marker='s',label='1%<σ<5%')     #1%<σ(H0)<5%
    if(second[0]>0):
        a2=plt.scatter(redshift1[[second]],M1[[second]],s=180,color='yellow',edgecolor='orange',marker='*',label='0.5%<σ<1%')                         #0.5%<σ(H0)<1%
    if(first[0]>0):
        a1=plt.scatter(redshift1[[first]],M1[[first]],s=80,color='r',marker='d',label='σ<0.5%')                  #σ(H0)<0.5%
    if(four[0]>0):
        a4=plt.scatter(redshift1[[four]],M1[[four]],s=80,color='w',edgecolor='b',label='σ>5%')
    #plt.grid(True)
    font={'size':18} 
    plt.legend(handles=[a1,a2,a3,a4,a0],labels=['$σ_{H_0}<0.5\%$','$0.5\%<σ_{H_0}<1\%$','$1\%<σ_{H_0}<5\%$','$σ_{H_0}>5\%$','$N_{host}>10^6$'],prop=font)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)    
    plt.title('PopIII in network',fontdict={'size':30})
    #plt.title('PopIII w/ lensing noise in Taiji')
    plt.xlim(0,20)
    plt.ylim(2.3,10.0)
    plt.ylabel('$\log_{10}{[M_c/M_{\odot}]}$', fontsize = 30)
    plt.xlabel('Redshift ($z$)', fontsize = 30)
    plt.savefig("fiveyears/network/lens/network_PopIII_lens%d.pdf"%datanum)
    #plt.savefig("threeyears/taiji/lens/PopIIIlens%d.png"%datanum, dpi=300)
    plt.show() 



##############Q3d
'''
for i in range(0,num):
    for j in range(0,N):
        if(dN_dtdzdM2[i][j]!=0):
            dN_dtdzdM2[i][j]=math.log(dN_dtdzdM2[i][j],10)
        else:
            dN_dtdzdM2[i][j]=-8    
'''

#network wirh lensing  in five years
exec('first{}=[-1]'.format(0));exec('second{}=[-1]'.format(0));exec('third{}=[849]'.format(0));exec('four{}=[850]'.format(0));
exec('first{}=[-1]'.format(1));exec('second{}=[-1]'.format(1));exec('third{}=[867]'.format(1));exec('four{}=[-1]'.format(1));
exec('first{}=[-1]'.format(2));exec('second{}=[-1]'.format(2));exec('third{}=[865]'.format(2));exec('four{}=[872,854,837]'.format(2));
exec('first{}=[-1]'.format(3));exec('second{}=[-1]'.format(3));exec('third{}=[-1]'.format(3));exec('four{}=[861]'.format(3));
exec('first{}=[-1]'.format(4));exec('second{}=[859,867]'.format(4));exec('third{}=[-1]'.format(4));exec('four{}=[865]'.format(4));
exec('first{}=[-1]'.format(5));exec('second{}=[-1]'.format(5));exec('third{}=[-1]'.format(5));exec('four{}=[848]'.format(5));
exec('first{}=[-1]'.format(6));exec('second{}=[-1]'.format(6));exec('third{}=[-1]'.format(6));exec('four{}=[841,849]'.format(6));
exec('first{}=[-1]'.format(7));exec('second{}=[861]'.format(7));exec('third{}=[867]'.format(7));exec('four{}=[841,839]'.format(7));
exec('first{}=[-1]'.format(8));exec('second{}=[-1]'.format(8));exec('third{}=[-1]'.format(8));exec('four{}=[850,858]'.format(8));
exec('first{}=[-1]'.format(9));exec('second{}=[-1]'.format(9));exec('third{}=[-1]'.format(9));exec('four{}=[-1]'.format(9));

datanum=0               
for datanum in range(4,5):
    #read the data                                                     #nine-dimensional matrix(M_c,eta,iota,d_L,alpha,delta,t_c,phi_c,psi)
    redshift1=np.loadtxt("fiveyears/redshift%d.txt"%datanum)
    M1=np.loadtxt("fiveyears/Mc%d.txt"%datanum)                                 #Msun
    i=0
    while(i<1456):
        M1[i]=math.log(M1[i],10)
        i=i+1
    sns.set_style("white")
    #plt.figure(figsize=(12,10))    
    plt.gcf().set_size_inches(10, 7)
    [X, Y] = np.meshgrid(z, log_M)
    list=np.array([-4.0,-2.0,-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.1,0.2])
    CQ3d1=plt.contourf(X, Y, dN_dtdzdM2, list,cmap=plt.cm.Greys,alpha=0.5,extend='min')
    CQ3d= plt.contour(X, Y, dN_dtdzdM2, list,colors=['k','k','k','k','k','k'])
    plt.clabel(CQ3d, inline=1, inline_spacing=0,fontsize=11,fmt='%1.1f')
    cbar=plt.colorbar(CQ3d1)
    cbar.set_label('$\log{[dN/(dtdzdlogM_c)]}$',fontsize=20)
    cbar.ax.tick_params(labelsize=15)
    #cbar.set_ticks([-8,-6,-4,-2,0,2])
    #cbar.set_ticklabels(['$0$','$10^{-6}$','$10^{-4}$','$10^{-2}$','$10^0$','$10^2$'])
    a0=plt.scatter(redshift1[835:875],M1[835:875],s=10,color='b')
    exec('first=first{}'.format(datanum))
    exec('second=second{}'.format(datanum))
    exec('third=third{}'.format(datanum))
    exec('four=four{}'.format(datanum))
    print(first,second,third)
    if(third[0]>0):
        a3=plt.scatter(redshift1[[third]],M1[[third]],s=100,color='lightgreen',marker='s',label='1%<σ<5%')     #1%<σ(H0)<5%
    if(second[0]>0):
        a2=plt.scatter(redshift1[[second]],M1[[second]],s=180,color='yellow',edgecolor='orange',marker='*',label='0.5%<σ<1%')                         #0.5%<σ(H0)<1%
    if(first[0]>0):
        a1=plt.scatter(redshift1[[first]],M1[[first]],s=80,color='r',marker='d',label='σ<0.5%')                  #σ(H0)<0.5%
    if(four[0]>0):
        a4=plt.scatter(redshift1[[four]],M1[[four]],s=80,color='w',edgecolor='b',label='σ>5%')
    #plt.grid(True)
    font={'size':18} 
    plt.legend(handles=[a1,a2,a3,a4,a0],labels=['$σ_{H_0}<0.5\%$','$0.5\%<σ_{H_0}<1\%$','$1\%<σ_{H_0}<5\%$','$σ_{H_0}>5\%$','$N_{host}>10^6$'],prop=font)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)      
    plt.title('Q3d in network',fontdict={"size":30})
    #plt.title('Q3d w/ lensing noise in Taiji')
    plt.xlim(0,20)
    plt.ylim(3.0,10.0)
    plt.ylabel('$\log_{10}{[M_c/M_{\odot}]}$', fontsize = 30)
    plt.xlabel('Redshift ($z$)', fontsize = 30)
    plt.savefig("fiveyears/network/lens/network_Q3d_lens%d.pdf"%datanum)
    #plt.savefig("threeyears/taiji/lens/Q3dlens%d.png"%datanum,dpi=300 )
    plt.show() 


##########lensing   in  one year
############Q3nod
'''
for i in range(0,num):
    for j in range(0,N):
        if(dN_dtdzdM3[i][j]!=0):
            dN_dtdzdM3[i][j]=math.log(dN_dtdzdM3[i][j],10)
        else:
            dN_dtdzdM3[i][j]=-8
'''




#network with lensing in  1year
exec('first{}=[-1]'.format(0));exec('second{}=[-1]'.format(0));exec('third{}=[227]'.format(0));exec('four{}=[-1]'.format(0));
exec('first{}=[-1]'.format(1));exec('second{}=[-1]'.format(1));exec('third{}=[194,208]'.format(1));exec('four{}=[-1]'.format(1));
exec('first{}=[-1]'.format(2));exec('second{}=[-1]'.format(2));exec('third{}=[-1]'.format(2));exec('four{}=[238,260,233]'.format(2));
exec('first{}=[-1]'.format(3));exec('second{}=[-1]'.format(3));exec('third{}=[281]'.format(3));exec('four{}=[223,245]'.format(3));
exec('first{}=[-1]'.format(4));exec('second{}=[-1]'.format(4));exec('third{}=[247,191,226,196]'.format(4));exec('four{}=[-1]'.format(4));
exec('first{}=[902]'.format(5));exec('second{}=[1024]'.format(5));exec('third{}=[1298,1293,1031,1442,992]'.format(5));exec('four{}=[929,883,1223,939,1050,1102,959]'.format(5));
exec('first{}=[-1]'.format(6));exec('second{}=[-1]'.format(6));exec('third{}=[1428,898,1453,1219]'.format(6));exec('four{}=[993,1370,1355,1393,1400,1345]'.format(6));
exec('first{}=[-1]'.format(7));exec('second{}=[-1]'.format(7));exec('third{}=[1271,1309,1257,1338,1006,1417,1196]'.format(7));exec('four{}=[987,883,993,1028,894,1292,1287,884,1439,1047,1245]'.format(7));
exec('first{}=[182]'.format(8));exec('second{}=[-1]'.format(8));exec('third{}=[-1]'.format(8));exec('four{}=[228]'.format(8));
exec('first{}=[-1]'.format(9));exec('second{}=[-1]'.format(9));exec('third{}=[1341,1141]'.format(9));exec('four{}=[1074,1354,1029]'.format(9));

datanum=0               
for datanum in [4,8]:
    #read the data                                                     #nine-dimensional matrix(M_c,eta,iota,d_L,alpha,delta,t_c,phi_c,psi)
    redshift1=np.loadtxt("data/redshift%d.txt"%datanum)
    M1=np.loadtxt("data/Mc%d.txt"%datanum)                                 #Msun
    i=0
    while(i<291):
        M1[i]=math.log(M1[i],10)
        i=i+1    
    sns.set_style("white")
    #plt.figure(figsize=(12,10))      
    plt.gcf().set_size_inches(10, 7)
    [X, Y] = np.meshgrid(z, log_M)
    list=np.array([-4,-3,-2.0,-1.0,0.0,0.5,1.0,])
    CQ3nod1=plt.contourf(X, Y, dN_dtdzdM3, list,cmap=plt.cm.Greys,alpha=0.5,extend='min')
    CQ3nod= plt.contour(X, Y, dN_dtdzdM3, list,colors=['k','k','k','k','k','k'])
    plt.clabel(CQ3nod, inline=1, inline_spacing=0,fontsize=8,fmt='%1.1f')

    cbar=plt.colorbar(CQ3nod1)
    cbar.set_label('$\log{[dN/(dtdzdlogM_c)]}$',fontsize=20)
    cbar.ax.tick_params(labelsize=15)
    #cbar.set_ticks([-8,-6,-4,-2,0,2])
    #cbar.set_ticklabels(['$0$','$10^{-6}$','$10^{-4}$','$10^{-2}$','$10^0$','$10^2$'])
    a0=plt.scatter(redshift1[175:291],M1[175:291],s=10,color='b')
    exec('first=first{}'.format(datanum))
    exec('second=second{}'.format(datanum))
    exec('third=third{}'.format(datanum))
    exec('four=four{}'.format(datanum))
    print(first,second,third)
    if(third[0]>0):
        a3=plt.scatter(redshift1[[third]],M1[[third]],s=100,color='lightgreen',marker='s',label='1%<σ<5%')     #1%<σ(H0)<5%
    if(second[0]>0):
        a2=plt.scatter(redshift1[[second]],M1[[second]],s=180,color='yellow',edgecolor='orange',marker='*',label='0.5%<σ<1%')                         #0.5%<σ(H0)<1%
    if(first[0]>0):
        a1=plt.scatter(redshift1[[first]],M1[[first]],s=80,color='r',marker='d',label='σ<0.5%')                  #σ(H0)<0.5%    
    if(four[0]>0):
        a4=plt.scatter(redshift1[[four]],M1[[four]],s=80,color='w',edgecolor='b',label='σ>5%')
    #plt.grid(True)
    font={'size':18}    
    plt.legend(handles=[a1,a2,a3,a4,a0],labels=['$σ_{H_0}<0.5\%$','$0.5\%<σ_{H_0}<1\%$','$1\%<σ_{H_0}<5\%$','$σ_{H_0}>5\%$','$N_{host}>10^6$'],prop=font)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)    
   
    #plt.title('Q3nod in network',fontdict={'size': 20})
    plt.title('Q3nod w/ lensing noise in network',fontdict={'size': 22})
    plt.xlim(0,20)
    plt.ylim(3.5,10.0)
    plt.ylabel('$\log_{10}{[M_c/M_{\odot}]}$', fontsize = 30)
    plt.xlabel('Redshift ($z$)', fontsize = 30)
    plt.savefig("data/network/lens/Q3nodlens%d_1yr_net.pdf"%datanum)
    #plt.savefig("fiveyears/taiji/lens/Q3nodlens%d.png"%datanum, dpi=300)
    plt.show() 

##########lensing   in  one year
############Q3nod
'''
for i in range(0,num):
    for j in range(0,N):
        if(dN_dtdzdM3[i][j]!=0):
            dN_dtdzdM3[i][j]=math.log(dN_dtdzdM3[i][j],10)
        else:
            dN_dtdzdM3[i][j]=-8
'''




#network without lensing in  1year
exec('first{}=[-1]'.format(0));exec('second{}=[-1]'.format(0));exec('third{}=[227]'.format(0));exec('four{}=[-1]'.format(0));
exec('first{}=[-1]'.format(1));exec('second{}=[-1]'.format(1));exec('third{}=[194,208]'.format(1));exec('four{}=[-1]'.format(1));
exec('first{}=[-1]'.format(2));exec('second{}=[-1]'.format(2));exec('third{}=[-1]'.format(2));exec('four{}=[238,260,233]'.format(2));
exec('first{}=[-1]'.format(3));exec('second{}=[-1]'.format(3));exec('third{}=[281]'.format(3));exec('four{}=[223,245]'.format(3));
exec('first{}=[247,191,196]'.format(4));exec('second{}=[222,226]'.format(4));exec('third{}=[179,212,185,207,210,261,225,218,220,229,254,242,211,224]'.format(4));exec('four{}=[284,186,264,241,215,187,180,275,176,237,199,184,217]'.format(4));
exec('first{}=[902]'.format(5));exec('second{}=[1024]'.format(5));exec('third{}=[1298,1293,1031,1442,992]'.format(5));exec('four{}=[929,883,1223,939,1050,1102,959]'.format(5));
exec('first{}=[-1]'.format(6));exec('second{}=[-1]'.format(6));exec('third{}=[1428,898,1453,1219]'.format(6));exec('four{}=[993,1370,1355,1393,1400,1345]'.format(6));
exec('first{}=[-1]'.format(7));exec('second{}=[-1]'.format(7));exec('third{}=[1271,1309,1257,1338,1006,1417,1196]'.format(7));exec('four{}=[987,883,993,1028,894,1292,1287,884,1439,1047,1245]'.format(7));
exec('first{}=[182]'.format(8));exec('second{}=[246]'.format(8));exec('third{}=[241,237,271,232,197,216,243,176,213,288,263,272]'.format(8));exec('four{}=[178,209,282,228,264,210,183]'.format(8));
exec('first{}=[-1]'.format(9));exec('second{}=[-1]'.format(9));exec('third{}=[1341,1141]'.format(9));exec('four{}=[1074,1354,1029]'.format(9));

datanum=0               
for datanum in [4,8]:
    #read the data                                                     #nine-dimensional matrix(M_c,eta,iota,d_L,alpha,delta,t_c,phi_c,psi)
    redshift1=np.loadtxt("data/redshift%d.txt"%datanum)
    M1=np.loadtxt("data/Mc%d.txt"%datanum)                                 #Msun
    i=0
    while(i<291):
        M1[i]=math.log(M1[i],10)
        i=i+1    
    sns.set_style("white")
    #plt.figure(figsize=(12,10))      
    plt.gcf().set_size_inches(10, 7)
    [X, Y] = np.meshgrid(z, log_M)
    list=np.array([-4,-3,-2.0,-1.0,0.0,0.5,1.0,])
    CQ3nod1=plt.contourf(X, Y, dN_dtdzdM3, list,cmap=plt.cm.Greys,alpha=0.5,extend='min')
    CQ3nod= plt.contour(X, Y, dN_dtdzdM3, list,colors=['k','k','k','k','k','k'])
    plt.clabel(CQ3nod, inline=1, inline_spacing=0,fontsize=8,fmt='%1.1f')

    cbar=plt.colorbar(CQ3nod1)
    cbar.set_label('$\log{[dN/(dtdzdlogM_c)]}$',fontsize=20)
    cbar.ax.tick_params(labelsize=15)
    #cbar.set_ticks([-8,-6,-4,-2,0,2])
    #cbar.set_ticklabels(['$0$','$10^{-6}$','$10^{-4}$','$10^{-2}$','$10^0$','$10^2$'])
    a0=plt.scatter(redshift1[175:291],M1[175:291],s=10,color='b')
    exec('first=first{}'.format(datanum))
    exec('second=second{}'.format(datanum))
    exec('third=third{}'.format(datanum))
    exec('four=four{}'.format(datanum))
    print(first,second,third)
    if(third[0]>0):
        a3=plt.scatter(redshift1[[third]],M1[[third]],s=100,color='lightgreen',marker='s',label='1%<σ<5%')     #1%<σ(H0)<5%
    if(second[0]>0):
        a2=plt.scatter(redshift1[[second]],M1[[second]],s=180,color='yellow',edgecolor='orange',marker='*',label='0.5%<σ<1%')                         #0.5%<σ(H0)<1%
    if(first[0]>0):
        a1=plt.scatter(redshift1[[first]],M1[[first]],s=80,color='r',marker='d',label='σ<0.5%')                  #σ(H0)<0.5%    
    if(four[0]>0):
        a4=plt.scatter(redshift1[[four]],M1[[four]],s=80,color='w',edgecolor='b',label='σ>5%')
    #plt.grid(True)
    font={'size':15}    
    plt.legend(handles=[a1,a2,a3,a4,a0],labels=['$σ_{H_0}<0.5\%$','$0.5\%<σ_{H_0}<1\%$','$1\%<σ_{H_0}<5\%$','$σ_{H_0}>5\%$','$N_{host}>10^6$'],prop=font)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)  
   
    #plt.title('Q3nod in network',fontdict={'size': 20})
    plt.title('Q3nod w/o lensing noise in network',fontdict={'size': 22})
    plt.xlim(0,20)
    plt.ylim(3.5,10.0)
    plt.ylabel('$\log_{10}{[M_c/M_{\odot}]}$', fontsize = 30)
    plt.xlabel('Redshift ($z$)', fontsize = 30)
    plt.savefig("data/network/nolens/Q3nod%d_1yr_net_nolens.pdf"%datanum)
    #plt.savefig("fiveyears/taiji/lens/Q3nodlens%d.png"%datanum, dpi=300)
    plt.show() 

