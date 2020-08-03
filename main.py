import numpy as np 
import matplotlib.pyplot as plt 
import math 
from mpl_toolkits.mplot3d import Axes3D

# parameters 
T = 1 
S_d = 30 
S_u = 83 
X_d = math.log(S_d) 
X_u = math.log(S_u) 
K = 50 
r = 0.08 
sigma = 0.25 
alpha = 1.5


# discretize 
N = 200 # time
M = 300 # spatial 
dt = T/N 
dx = (X_u-X_d)/M 
t = np.linspace(0,T,N+1) 
X = np.linspace(X_d,X_u,M+1) 
S = np.exp(X)

# midumn 
v = -sigma**alpha/math.cos(alpha*math.pi/2)/2 
zeta = v*dt/(2*dx**alpha) 
xi = dt*(r-v)/(4*dx) 
eta = dt*r/2 

# Matrix 
g = np.zeros((M+1,1)) 
omega = np.zeros((M+1,1)) 
g[0,0] = 1 
omega[0,0] = alpha*g[0,0]/2 
for k in range(1,M+1): 
    g[k,0] = (1-(alpha+1)/k)*g[k-1,0] 
    omega[k,0] = alpha*g[k,0]/2+(2-alpha)*g[k-1,0]/2 

A = np.zeros((M-1,M-1)) 
for i in range(M-1): 
    for j in range(M-1): 
        if i == j: 
            A[i,j] = omega[1,0]
        elif i-j ==1:
            A[i,j] = omega[2,0]
        elif i-j == -1:
            A[i,j] = omega[0,0]
        elif i-j>=2:
            A[i,j] = omega[i-j+1,0]
        else:
            A[i,j] = 0
B = np.diag(-1*np.ones((M-2,)),k=-1)+np.diag(np.ones((M-2,)),k=1) 

L = (1+eta)*np.eye(M-1)-zeta*A-xi*B 
R = (1-eta)*np.eye(M-1)+zeta*A+xi*B 

# calculate 
V = np.zeros((N+1,M+1)) 
for i in range(M+1): 
    if X_d+i*dx <= math.log(83):
        V[0,i] = max(math.exp(X_d+i*dx)-K,0)
    else:
        V[0,i] = 0
V[:,0:1] = np.zeros((N+1,1)) 
V[:,M:M+1] = np.zeros((N+1,1)) 
for i in range(1,N+1):
    V[i:i+1,1:M] = np.transpose(np.linalg.inv(L).dot(R).dot(V[i-1:i,1:M].T)) 

# picture
S, t = np.meshgrid(S, t)
fig1 = plt.figure(1)
ax = Axes3D(fig1)
ax.plot_surface(S, T-t, V, rstride=1, cstride=1, cmap='rainbow')
plt.xlabel('Price of stock')
plt.ylabel('Time')
ax.set_zlabel('Price of option')

#fig2 = plt.figure(2)
#plt.plot(np.exp(X), V[51,:]) 
#plt.plot(np.exp(X), V[int(N/4),:]) 
#plt.plot(np.exp(X), V[int(N/4),:]) 
#plt.plot(np.exp(X), V[int(N/6),:]) 
#plt.plot(np.exp(X), V[int(N/12),:]) 
#plt.plot(np.exp(X), V[int(N/26),:]) 
#plt.plot(np.exp(X), V[int(N/52),:]) 

plt.show() 