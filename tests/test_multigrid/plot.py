import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys

plt.ion()
j = int(sys.argv[1])
n = 1 + 2**(j+1)
ans = np.fromfile('sols/solution{:d}.dat'.format(j),dtype=np.float64).reshape(n,n)
res = np.fromfile('result.dat').reshape(n,n)
src = np.fromfile('source.dat').reshape(n,n)
grid = np.fromfile('grid.dat')
x = grid[:n]
y = grid[-n:]
xx,yy = np.meshgrid(x,y,indexing='ij')


ans += 1e-20

err = abs((res-ans))


print("Min / Max error {:.2e}, {:.2e}".format(abs(res-ans).flatten().min(), abs(res-ans).flatten().max()))
print("L2 error: {:.2e}".format(np.sqrt(sum(err.flatten()**2))))


fig,axes = plt.subplots(2,2,figsize=(15,8),sharex=True,sharey=True)

imgs = [None]*4
imgs[0]=axes[0,0].pcolormesh(xx,yy,ans)
imgs[1]=axes[0,1].pcolormesh(xx,yy,src)
imgs[2]=axes[1,0].pcolormesh(xx,yy,res)
imgs[3]=axes[1,1].pcolormesh(xx,yy,np.log10(err+1e-8))

for ax,im in zip(axes.flatten(),imgs):
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im, cax=cax, orientation='vertical')
    ax.set_aspect('equal')
axes[0,0].set_title('Answer')
axes[1,0].set_title('Result')
axes[1,1].set_title('Resid')
axes[0,1].set_title('Source')


fig.tight_layout()

plt.show()
plt.pause(2**31-1)
