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




##sensitivity noise curve of lensing in redshift of 0.3-3.0
Label=['z=0.3','z=0.6','z=0.9','z=1.2','z=1.5','z=1.8','z=2.1','z=2.4','z=2.7','z=3.0']
Label1=['z=0.3','z=0.6','z=0.9','z=1.2','z=1.5','z=1.8','z=2.1','z=2.4','z=2.7','z=3.0']
value=range(10)
colors=plt.get_cmap('jet')
cNorm=mpl.colors.Normalize(vmin=0,vmax=value[-1])
scalarMap=cm.ScalarMappable(norm=cNorm,cmap=colors)
#print(value,value[-1])
sns.set_style("white")
#plt.figure(figsize=(12,10))
plt.gcf().set_size_inches(7,5)


plt.plot([0,0],[0,0],linestyle=' ',label="Lensing noise")
for index in range(0,10):
    f_lens=np.loadtxt("lens/f_lens%d.txt"%index)         #  frequency of lensing noise
    n_lens=np.loadtxt("lens/lens%d.txt"%index)       #  noise of lensing  in redshift of 0.3-3.0
    for i in range(len(f_lens)):
       f_lens[i]=math.log(f_lens[i],10) 
       n_lens[i]=math.log(n_lens[i],10)
    colorVal=scalarMap.to_rgba(index)
    plt.plot(f_lens,n_lens,c=colorVal,linestyle='-',label=Label[index],linewidth=1)



#sensitivity noise curve of LISA
f_lisa=np.loadtxt("detector/f_LISA.txt")
s_lisa=np.loadtxt("detector/sen_LISA.txt")     

for i in range(len(f_lisa)):
    f_lisa[i]=math.log(f_lisa[i],10)
    s_lisa[i]=math.log(s_lisa[i],10)

plt.plot(f_lisa,s_lisa,color='k',linewidth='2',linestyle='-',label="LISA")



plt.plot([0,0],[0,0],linestyle=' ',label="signal")
#signal curve 
for index in range(0,10):
    f_sig=np.loadtxt("sig/f_sig%d.txt"%index)           # signal
    s_sig=np.loadtxt("sig/sig%d.txt"%index)
    #print(np.max(n_lens),np.min(n_lens))
    for i in range(len(f_sig)):
       f_sig[i]=math.log(f_sig[i],10)
       s_sig[i]=math.log(s_sig[i],10)
    colorVal=scalarMap.to_rgba(index)
    plt.plot(f_sig,s_sig,c=colorVal,linestyle='-',label=Label1[index],linewidth=3)






#sensitivity noise curve of taiji    
f_taiji=np.loadtxt("detector/f_taiji.txt")
s_taiji=np.loadtxt("detector/sen_Taiji.txt")

for i in range(len(f_taiji)):
    f_taiji[i]=math.log(f_taiji[i],10)
    s_taiji[i]=math.log(s_taiji[i],10)

plt.plot(f_taiji,s_taiji,color='k',linewidth='2',linestyle='--',label="Taiji")  




plt.xticks([-5,-4,-3,-2,-1,0],['$10^{-5}$','$10^{-4}$','$10^{-3}$','$10^{-2}$','$10^{-1}$','$10^{0}$'],fontsize=12)
plt.yticks([-21,-20,-19,-18,-17,-16,-15,-14,-13,-12],
           ['$10^{-21}$','$10^{-20}$','$10^{-19}$','$10^{-18}$','$10^{-17}$','$10^{-16}$','$10^{-15}$','$10^{-14}$','$10^{-13}$','$10^{-12}$'],
           fontsize=12)
plt.xlim(-5,0)
plt.ylim(-21,-12)
plt.grid(linestyle='--')
font={'size':9}
plt.legend(ncol=2,prop=font)
#plt.title('lensing noise')
plt.xlabel('Frequency[Hz]',fontsize=15)
plt.ylabel('$Strain[Hz^{-1/2}$]',fontsize=15)
plt.savefig("sensitivity_noise_curve.pdf")
plt.show()





#signal curve 
for index in range(0,2):
    f_sig=np.loadtxt("sig/f_sig%d.txt"%index)           # signal
    s_sig=np.loadtxt("sig/sig%d.txt"%index)
    #print(np.max(n_lens),np.min(n_lens))
    for i in range(len(f_sig)):
       f_sig[i]=math.log(f_sig[i],10)
       s_sig[i]=math.log(s_sig[i],10)
    colorVal=scalarMap.to_rgba(index)
    plt.plot(f_sig,s_sig,c=colorVal,linestyle='--',label=Label1[index],linewidth=2.5)


plt.xticks([-5,-4,-3,-2,-1,0],['$10^{-5}$','$10^{-4}$','$10^{-3}$','$10^{-2}$','$10^{-1}$','$10^{0}$'])
plt.yticks([-21,-20,-19,-18,-17,-16,-15,-14,-13,-12],
           ['$10^{-21}$','$10^{-20}$','$10^{-19}$','$10^{-18}$','$10^{-17}$','$10^{-16}$','$10^{-15}$','$10^{-14}$','$10^{-13}$','$10^{-12}$'])
plt.xlim(-5,0)
plt.ylim(-21,-12)
font={'size':8}
plt.grid(linestyle='--')
plt.legend(ncol=2)
#plt.title('lensing noise')
plt.xlabel('Frequency[Hz]')
plt.ylabel('$Strain[Hz^{-1/2}$]')
plt.savefig("signal.png",dpi=500)
plt.show()

