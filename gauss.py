import numpy as np
import ase
from ase.io import read 
import matplotlib.pyplot as plt
import sys

data=ase.io.read('pr_npt.xyz',':')
_,_,_,_,_,_,a,b,c=np.loadtxt('pr_npt.out',skiprows=1,unpack=True)

def get_bpa_coms(data,N=1):
    coms=[]
    for i in range(N):
        a=1040+(33*i) ; b=1040+(33*(i+1))
        bpa=data[a:b]
        coms.append(bpa.get_center_of_mass())
    return coms

def wrap(x,xc,L):
    dx=float(x)-float(xc)
    if dx > L/2:
        dx-=L
    elif dx < -L/2:
        dx+=L
    return dx

def add_gaussians(cnt,coms,N=1,sigma=1,height=1,cnt_or_graphene='cnt'):
    gaussians=np.zeros(len(cnt))
    Lx=cnt.get_cell()[0][0]
    Lz=cnt.get_cell()[-1][-1]
    for j in range(len(coms)):
        xc,yc,zc=coms[j]
        for i in range(len(cnt)):
            x,y,z=cnt.get_positions()[i]
            if cnt_or_graphene="cnt":
                dx,dy,dz=x-xc,y-yc,wrap(z,zc,L)
            elif cnt_or_graphene="graphene":
                dx,dy,dz=wrap(x,xc,Lx),y-yc,wrap(z,zc,L)
            gaussians[i]+=height*np.e**(-(((dx**2)+(dy**2)+(dz**2))/(2*(sigma**2))))
    return gaussians

def frame_gaussians(data,a,b,c,Nframe=1,N=1,sigma=1,height=1,cnt_or_graphene='cnt'):
    new_gaussians=[]
    for i in range(Nframe):
        current=data[i]
        current.set_cell([a[i],b[i],c[i]])
        cnt=current[:1040]
        coms=get_bpa_coms(current,N)
        new_gaussians.append(add_gaussians(cnt,coms,N,sigma,height,cnt_or_graphene))
    return new_gaussians

def get_sums(gauss):
    sums=[]
    for i in gauss:
        sums.append(np.sum(i))
    return sums

def block_avg(sums,n_block=100):
    blocksize=int(np.floor(len(sums)/n_block))
    b_avg=[]

    for i in range(n_block):
        i1=blocksize*i ; i2=blocksize*(i+1)
        block=sums[i1:i2]
        b_avg.append(np.average(block))
        if i==n_block-1:
            i1=blocksize*i
            block=sums[i1:]
            b_avg.append(np.average(block))
    avg=np.average(b_avg)
    err=np.sqrt(np.sum((b_avg-avg)**2)/n_block)
    return avg,err

gauss=frame_gaussians(data,a,b,c,Nframe=len(data),N=int(sys.argv[-1]),sigma=3,height=1)
sums=get_sums(gauss)

np.savetxt('coverages.out',np.array(sums))

avg,err=block_avg(sums)
print(avg,err)
