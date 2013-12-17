import numpy as np
import matplotlib.pyplot as plt
import sys, os

# Kernel height
h   = 1.
# Maximum interaction distance (Kernel dependant)
sep = 2.
# Number of points considered
N   = 10000
# x coordinates
x   = np.linspace(-sep,sep,N)
# Evaluate the kernel and the gradient
w   = np.zeros(N)
dw  = np.zeros(N)
for i in range(0,N):
	w[i]  =  7./(64. * np.pi * h**2.) * (1.+2.*np.abs(x[i])) * (2.-np.abs(x[i]))**4.
	dw[i] = 70./(64. * np.pi * h**3.) * x[i] * (2.-np.abs(x[i]))**3.
# Plot results
plt.xlabel(r'$\vert x \vert / h$', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
wLine,  = plt.plot(x, w,  label=r'$W(\vert x \vert; h)$', lw=2)
dwLine, = plt.plot(x, dw, label=r'$\nabla W(\vert x \vert; h)$', lw=2)
plt.legend(loc='best', prop={'size':20})
plt.grid()
pathname = os.path.dirname(sys.argv[0])
plt.savefig(os.path.join(pathname, 'wendland2D.pdf'), format='pdf')