import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation



def Param():
    
    L= 20.
    dx = 0.01
    n = int(L/dx)
    X = np.linspace(0,L,n)
    x = X/L
    dt = 10**(-5)
    
    tmax = 2.
    t = int( tmax / dt )
    
    hbar = 1.
    m = 1.
    xs = 5.
    kx = 20.
    sigma = 1.
    
    E = hbar**2 * kx**2 / (2 * m)
    return L, dx, X, x, xs, kx, n, dt, hbar, m, sigma, t, tmax, E

    
def Phi(sigma, X, dx, xs, kx, n, L):
    PhiR = np.exp(-((X-xs)**2 ) / ( 2. * sigma**2 ) ) * np.cos( kx * X )
    PhiI = np.exp(-((X-xs)**2 ) / ( 2. * sigma**2 ) ) * np.sin( kx * X )
    C = Normal(PhiR, PhiI, sigma, X, dx, xs, kx, n, L)
    PhiR = C * PhiR
    PhiI = C * PhiI
    PhiI[0] = PhiI[n-1] = PhiR[0] = PhiR[n-1] = 0
    return PhiR, PhiI


def Normal(PhiR, PhiI, sigma, X, dx, xs, kx, n, L):
    Prob = PhiR**2 + PhiI**2
    c = 0
    for i in range(len(X)):
        c = c + Prob[i] * dx
    C = 1/c
    return C

def Phit(PhiR, PhiI, V,  dx, dt, n, hbar, m):
    PhiI[1:n-1] = PhiI[1:n-1] - dt * ( PhiR[1:n-1] * V[1:n-1] / hbar - hbar * (PhiR[2:n] - 2 * PhiR[1:n-1] + PhiR[0:n-2]) / ( 2 * m * dx**2 ))
    PhiR[1:n-1] = PhiR[1:n-1] + dt * ( PhiI[1:n-1] * V[1:n-1] / hbar - hbar * (PhiI[2:n] - 2 * PhiI[1:n-1] + PhiI[0:n-2]) / ( 2 * m * dx**2 ))
    return PhiI, PhiR

def Pot_carre(E, n):
    V = [0]*n
    v = [0]*n
    z=int(n/2)
    for i in range(n):
        if i < z or i>z+500:
            V[i] = 0
            v[i] = V[i] / E
        else:
            V[i] = 250
            v[i] = V[i] / E
    return V, v

def Pot_tunnel(E, n):
    V = [0]*n
    v = [0]*n
    z=int(n/2)
    for i in range(n):
        if i < z or i > z+15:
            V[i] = 0
            v[i] = V[i] / E
        else:
            V[i] = 250
            v[i] = V[i] / E
    return V, v

def Pot_period(E, n):
    V = [0]*n
    v = [0]*n
    z=int(n/2)
    for i in range(n):
        if  z < i < z+15 or 11/10 * z < i < 11/10 * z+15 or  12/10 * z < i <  12/10 * z+15 or 13/10 * z < i <  13/10 * z+15:
            V[i] = 250
            v[i] = V[i] / E
        else:
            V[i] = 0
            v[i] = V[i] / E
    return V, v

def Pot_nul(E, n):
    V = [0]*n
    v = [0]*n
    return V, v


def Prob(PhiR, PhiI, V,  dx, dt, t, n, hbar, m):
    ProbTot = [PhiR**2 + PhiI**2]
    for j in range(t):
        PhiI, PhiR = Phit(PhiR, PhiI, V, dx, dt, n, hbar, m)
        ProbTot.append(PhiR**2 + PhiI**2)
    return ProbTot
    
    
def animat(Pot):
    
    L, dx, X, x, xs, kx, n, dt, hbar, m, sigma, t, tmax, E = Param()
    
    PhiR, PhiI = Phi(sigma, X, dx, xs, kx, n, L)
    
    V, v = Pot(E, n)
    
    ProbTot = Prob(PhiR, PhiI, V,  dx, dt, t, n, hbar, m)
    ProbTot = ProbTot[0:-1:100]
    
    fig = plt.figure()
    plt.xlabel('x')
    plt.plot(X,v)
    line, = plt.plot([],[]) 
    plt.xlim(0, L)
    plt.ylim(0,1)
    
    def init():
        line.set_data([],[])
        return line,
    
    
    def anim(i, ProbTot):
        ProbTot = ProbTot[i]
        line.set_data(X, ProbTot)
        return line,
    
    l = len(ProbTot)
    ani = animation.FuncAnimation(fig, anim, init_func=init, frames= l, fargs=(ProbTot, ), blit=True, interval=10, repeat=True)
    
    return ani

def enreg_ani(Pot):
    
    animation.rcParams['animation.ffmpeg_path'] = r'C:\ffmpeg\bin\ffmpeg.exe'
    
    mywriter = animation.FFMpegFileWriter(fps=25,codec="libx264")
    
    ani = animat(Pot)
    
    ani.save("Pot_carre_tmax=1.mp4", writer=mywriter)
    
    return ani
    