# coding: utf-8

# In[1]:

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors


# In[2]:

#skw = np.loadtxt("S_k_w.txt",usecols=(2,4,5))

skw1 = np.loadtxt("Akw.txt",usecols=(0,5,6))
skw2 = np.loadtxt("Akw.txt",usecols=(0,5,7))
skw3 = np.loadtxt("Akw.txt",usecols=(0,5,8))
skw4 = np.loadtxt("Akw.txt",usecols=(0,5,9))
skw5 = np.loadtxt("Akw.txt",usecols=(0,5,10))
skw6 = np.loadtxt("Akw.txt",usecols=(0,5,11))

#skw_av = np.loadtxt("Skw_OnAveragedConfs.txt",usecols=(2,4,6))

#print(skw1)
#print(skw2)
#print(skw3)

# In[4]:

x = skw1[:,0]
y = (skw1[:,1])

#z=100*( (1.0/20.0)*( ( 1.0*skw1[:,2])  + (1.0*skw2[:,2]) + (1.0*skw3[:,2])  + (1.0*skw4[:,2])  + ( 1.0*skw5[:,2])  + (1.0*skw6[:,2]) + (1.0*skw7[:,2])  + (1.0*skw8[:,2]) + ( 1.0*skw9[:,2])  + (1.0*skw10[:,2]) + (1.0*skw11[:,2])  + (1.0*skw12[:,2]) + ( 1.0*skw13[:,2])  + (1.0*skw14[:,2]) + (1.0*skw15[:,2])  + (1.0*skw16[:,2]) + ( 1.0*skw17[:,2])  + (1.0*skw18[:,2]) + (1.0*skw19[:,2]) + (1.0*skw20[:,2])  )  - 1.0*skw_av[:,2] )

z= skw1[:,2] + skw2[:,2] + skw3[:,2] + skw4[:,2] + skw5[:,2] + skw6[:,2]
#z= skw2[:,2] + skw5[:,2] 
#z= skw1[:,2] + skw2[:,2] + skw3[:,2]
ly = np.int(y.max())
lx = np.int(x.max())

print(np.size(x))
ssp = np.zeros((ly+1,lx+1))
for k in xrange(np.size(x)):
    a = np.int(x[k])
    b= np.int(y[k])
    ssp[b,a] = (z[k])
z_min, z_max = z.min(), z.max()
#print(z_min,z_max)

ssp_new = np.zeros((ly+1,lx+1))
for k in xrange(np.size(x)):

        a = np.int(x[k])
        b= np.int(y[k])

fig = plt.figure()
#fig, ax = plt.subplots(2, 1)
plt.ylabel('$\omega$',size=28, labelpad=0)
plt.xlabel('$q$',size=28, labelpad=0)
xticks_labels=['$0$', '$\pi$', '$2\pi$']
plt.xticks([0.5, 8.5, 16],xticks_labels, size=18)
yticks_labels=['']

plt.yticks([50],yticks_labels, size=24)
plt.xlim(xmin=0,xmax=26)
plt.ylim(ymin=0,ymax=199)
plt.pcolormesh(ssp,cmap = plt.cm.get_cmap('bwr') , norm = colors.Normalize(vmin=0,vmax=200))  
#plt.pcolormesh(ssp,cmap = plt.cm.get_cmap('wbr') , norm=colors.LogNorm(vmin=0.001, vmax=10))

#pcm = ax[1].pcolor(ssp, cmap='PuBu_r')
#fig.colorbar(extend='max')
plt.colorbar()
plt.show()

#fig.colorbar(pcm, ax=ax[0], extend='max')
fig.savefig('trial')

# In[ ]

