import numpy as np
import matplotlib as mpl
import math
from scipy import integrate
from scipy.optimize import fsolve,root
import sys
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import norm
import sympy
from time import * 
from scipy.spatial import ConvexHull
from tqdm import tqdm
from scipy import interpolate
import seaborn as sns
#plt.switch_backend('agg')

'''
#read the data
f_lens=np.loadtxt("f_0.txt")        # frequency of lensing noise
f_det=np.loadtxt("f_1.txt")
n_lens0=np.loadtxt("lens_0.txt")       #  noise of lensing  in redshift of 0.3
n_lens1=np.loadtxt("lens_1.txt")        #0.6
n_lens2=np.loadtxt("lens_2.txt")           #0.9
n_lens3=np.loadtxt("lens_3.txt")        #1.2
n_lens4=np.loadtxt("lens_4.txt")       #1.5
n_lens5=np.loadtxt("lens_5.txt")      #1.8
n_lens6=np.loadtxt("lens_6.txt")      #2.1
n_lens7=np.loadtxt("lens_7.txt")      #2.4
n_lens8=np.loadtxt("lens_8.txt")     #2.7
n_lens9=np.loadtxt("lens_9.txt")      #3
s_lisa=np.loadtxt("S_LISA.txt")
s_taiji=np.loadtxt("S_Taiji.txt")

###noise spectral
print(np.max(f_lens),np.min(f_lens))
print(np.max(f_det),np.min(f_det))
print(np.max(s_lisa),np.min(s_lisa))
print(np.max(s_taiji),np.min(s_taiji))
for i in range(len(f_lens)):
    f_lens[i]=math.log(f_lens[i],10)
for i in range(len(f_det)):
    f_det[i]=math.log(f_det[i],10)
    s_lisa[i]=math.log(s_lisa[i],10)
    s_taiji[i]=math.log(s_taiji[i],10)




colors=['magenta','brown','red','orange','yellow','green','cyan',
      'blue','purple','pink']
Label=['z=0.3','z=0.6','z=0.9','z=1.2','z=1.5','z=1.8','z=2.1','z=2.4','z=2.7','z=3.0']
plt.gcf().set_size_inches(7,5)
for index in range(0,10):
    n_lens=np.loadtxt("lens_%d.txt"%index)       #  noise of lensing  in redshift of 0.3-3.0
    #print(np.max(n_lens),np.min(n_lens))
    for i in range(len(n_lens)):
        n_lens[i]=math.log(n_lens[i],10)
    plt.plot(f_lens,n_lens,c=colors[index],label=Label[index])

plt.plot(f_det,s_lisa,color='k',linewidth='1',linestyle='-',label="LISA")
plt.plot(f_det,s_taiji,color='k',linewidth='1',linestyle='--',label="Taiji")    

plt.xticks([-5,-4,-3,-2,-1,0],['$10^{-5}$','$10^{-4}$','$10^{-3}$','$10^{-2}$','$10^{-1}$','$10^{0}$'])
plt.yticks([-21,-20,-19,-18,-17,-16,-15,-14],
           ['$10^{-21}$','$10^{-20}$','$10^{-19}$','$10^{-18}$','$10^{-17}$','$10^{-16}$','$10^{-15}$','$10^{-14}$'])
plt.grid(linestyle='--')
plt.legend()
#plt.title('lensing noise')
plt.xlabel('Frequency[Hz]')
plt.ylabel('$Spectral\ density[Hz^{-1/2}$]')
plt.savefig("noise.png",dpi=500)
plt.show()


value=range(10)
colors=plt.get_cmap('jet')
cNorm=mpl.colors.Normalize(vmin=0,vmax=value[-1])
scalarMap=cm.ScalarMappable(norm=cNorm,cmap=colors)
#print(value,value[-1])
plt.gcf().set_size_inches(7,5)
for index in range(0,10):
    n_lens=np.loadtxt("lens_%d.txt"%index)       #  noise of lensing  in redshift of 0.3-3.0
    #print(np.max(n_lens),np.min(n_lens))
    for i in range(len(n_lens)):
        n_lens[i]=math.log(n_lens[i],10)
    colorVal=scalarMap.to_rgba(index)
    plt.plot(f_lens,n_lens,c=colorVal,label=Label[index],linewidth=0.5)

plt.plot(f_det,s_lisa,color='k',linewidth='1',linestyle='-',label="LISA")
plt.plot(f_det,s_taiji,color='k',linewidth='1',linestyle='--',label="Taiji")    

plt.xticks([-5,-4,-3,-2,-1,0],['$10^{-5}$','$10^{-4}$','$10^{-3}$','$10^{-2}$','$10^{-1}$','$10^{0}$'])
plt.yticks([-21,-20,-19,-18,-17,-16,-15,-14],
           ['$10^{-21}$','$10^{-20}$','$10^{-19}$','$10^{-18}$','$10^{-17}$','$10^{-16}$','$10^{-15}$','$10^{-14}$'])

plt.grid(linestyle='--')
plt.legend()
#plt.title('lensing noise')
plt.xlabel('Frequency[Hz]')
plt.ylabel('$Spectral\ density[Hz^{-1/2}$]')
plt.savefig("noise1.png",dpi=500)
plt.show()

'''

####network with lensing
N=17
y1=np.linspace(0,N,N)         #number of events
x1=np.array([63.48348348,67.4734947,66.57331466,60.50251256,67.85357071,
             67.76835367,63.2160804,67.5675675,
             65.70570571,66.28140704,67.61130565,67.73334667,66.72672673,58.89447236,66.84684685,66.36636637,63.903903])                 
#142,222,446,532,590,
#859,865,867,
#971,992,1016,1181,1197,1203,1287,1408,1445
xerror=np.array([[23.48348348,6.246246246],
                [2.032406481,2.032406481],
                [6.913382677,6.565313063],
                [20.50251256,8.442211055],
                [0.280056011,0.26005201],
                [0.369673935,0.364872975],
                [23.2160804,6.030150754],
                [0.660660661,0.640640641],
                [5.405405405,4.204204204],
                [4.020100503,3.417085427],
                [0.624312156,0.591295648],      
                [0.075015003,0.075315063],
                [3.243243243,2.822822823],
                [18.89447236,9.648241206],
                [2.702702703,2.282282282],
                [3.543543544,2.822822823],
                [23.9039039,6.126126126]])
xerror=np.transpose(xerror)
#print(np.array(xerror[:,0]))
'''
fig,ax1=plt.subplots()
ax2=ax1.twinx()
ax2.spines["left"].set_position("axes",1.2)
'''

sns.set_style("white")
#plt.figure(figsize=(12,10))
value=range(N)
colors=plt.get_cmap('jet')
cNorm=mpl.colors.Normalize(vmin=0,vmax=value[-1])
scalarMap=cm.ScalarMappable(norm=cNorm,cmap=colors)
plt.gcf().set_size_inches(9,5)
for i in range(0,N):
    colorVal=scalarMap.to_rgba(i)
    plt.scatter(x1[i],y1[i],marker='o',color=colorVal,s=10)
    plt.errorbar(x1[i],y1[i],xerr=[[xerror[0][i]],[xerror[1][i]]],c=colorVal)



plt.yticks(y1,['PopIII-142','PopIII-222','PopIII-446','PopIII-532','PopIII-590','Q3d-859','Q3d-865','Q3d-867',
               'Q3nod-971','Q3nod-992','Q3nod-1016','Q3nod-1181','Q3nod-1197','Q3nod-1203','Q3nod-1287','Q3nod-1408','Q3nod-1445'],fontsize=12)
plt.xticks(fontsize=20)
plt.xlabel('$H_0$[$km\ s^{-1}$$Mpc^{-1}$]',fontsize=25)
plt.title('Network',fontsize=30)
plt.axvline(x=67.4,linestyle='--',color='grey')
plt.axvline(x=74.03,linestyle='--',color='orange')
plt.gca().add_patch(plt.Rectangle((67.4-0.5,-1), 0.5+0.5, N+2,facecolor="grey",alpha=0.5,label='CMB'))   #CMB: 67.4+-0.5
plt.gca().add_patch(plt.Rectangle((74.03-1.42,-1), 1.42+1.42, N+2,facecolor="lightblue",alpha=0.5,label='SNIa(SH0ES)')) 
plt.xlim(35,75)
font = { 'size': 15}
plt.legend(prop=font)
plt.savefig("fig/net_lens.pdf")
plt.show()



#### taiji with lensing
N=7
y2=np.linspace(0,N,N)
x2=np.array([67.64152831,67.80780781,67.44772386,67.80780781,67.57778889,67.73754751,61.30653266])      
#222,590,
#859,867,
#1016,1181,1197
xerror1=np.array([[5.833166633,4.656931386],
                  [1.021021021,1.081081081],
                  [0.980490245,0.825412706],
                  [0.900900901,0.900900901],
                  [0.84042021,0.745372686],
                  [0.104020804,0.104020804],
                  [21.30653266,8.442211055]])
xerror1=np.transpose(xerror1)
#print(xerror1,xerror1[0][1])

sns.set_style("white")
#plt.figure(figsize=(12,10))
value=range(N)
c=plt.get_cmap('jet')
cNorm=mpl.colors.Normalize(vmin=0,vmax=value[-1])
scalarMap=cm.ScalarMappable(norm=cNorm,cmap=c)
plt.gcf().set_size_inches(9,5)
for i in range(0,N):
    colorVal=scalarMap.to_rgba(i)
    plt.scatter(x2[i],y2[i],marker='o',color=colorVal,s=10)
    plt.errorbar(x2[i],y2[i],xerr=[[xerror1[0][i]],[xerror1[1][i]]],c=colorVal)



plt.yticks(y2,['PopIII-22','PopIII-590','Q3d-859','Q3d-867','Q3nod-1016','Q3nod-1181','Q3nod-1197'],fontsize=12)
plt.xticks(fontsize=20)
plt.xlabel('$H_0$[$km\ s^{-1}$$Mpc^{-1}$]',fontsize=25)
plt.title('Taiji',fontsize=30)
plt.axvline(x=67.4,linestyle='--',color='grey')
plt.axvline(x=74.03,linestyle='--',color='orange')
plt.gca().add_patch(plt.Rectangle((67.4-0.5,-1), 0.5+0.5, N+1.5,facecolor="grey",alpha=0.5,label='CMB'))   #CMB: 67.4+-0.5
plt.gca().add_patch(plt.Rectangle((74.03-1.42,-1), 1.42+1.42, N+1.5,facecolor="lightblue",alpha=0.5,label='SNIa(SH0ES)')) 
plt.xlim(35,80)
font={'size':15}
plt.legend(loc=0,prop=font)
plt.savefig("fig/taiji_lens.pdf")
plt.show()




### the best in network or taiji with lensing 
N=6
y3=np.linspace(0,N,N)
x3=np.array([[67.85357071],[67.76835367],[67.56756757],[67.61130565],[67.73334667],
             [67.73754751]])                                
#network590,859,867,1016,1181    taiji1181
xerror=np.array([[0.280056011,0.26005201],
                 [0.369673935,0.364872975],
                 [0.660660661,0.640640641],
                 [0.624312156,0.591295648],
                 [0.075015003,0.075315063],
                 [0.104020804,0.104020804]])
xerror=np.transpose(xerror)

sns.set_style("white")
#plt.figure(figsize=(12,10))
value=range(N)
c=plt.get_cmap('jet')
cNorm=mpl.colors.Normalize(vmin=0,vmax=value[-1])
scalarMap=cm.ScalarMappable(norm=cNorm,cmap=c)
plt.gcf().set_size_inches(15,8)
for i in range(0,N):
    colorVal=scalarMap.to_rgba(i)
    plt.scatter(x3[i],y3[i],marker='o',color=colorVal,s=10)
    plt.errorbar(x3[i],y3[i],xerr=[[xerror[0][i]],[xerror[1][i]]],c=colorVal)



plt.yticks(y3,['Network+PopIII-590','Network+Q3d-859','Network+Q3d-867','Network+Q3nod-1016','Network+Q3nod-1181','Taiji+Q3nod-1181'],
           fontsize=13,rotation=10)
plt.xticks(fontsize=30)
#plt.tick_params(labelsize=9)
plt.xlabel('$H_0$[$km\ s^{-1}$$Mpc^{-1}$]',fontsize=30)

#plt.title('taiji with lensing noise')
plt.axvline(x=67.4,linestyle='--',color='grey')
plt.axvline(x=74.03,linestyle='--',color='orange')
plt.gca().add_patch(plt.Rectangle((67.4-0.5,-1), 0.5+0.5, N+1.5,facecolor="grey",alpha=0.5,label='CMB'))   #CMB: 67.4+-0.5
plt.gca().add_patch(plt.Rectangle((74.03-1.42,-1), 1.42+1.42, N+1.5,facecolor="lightblue",alpha=0.5,label='SNIa(SH0ES)')) 
plt.xlim(65,76)
font={'size':13}
plt.legend(loc=2,prop=font)
plt.savefig("fig/best.pdf")
plt.show()











