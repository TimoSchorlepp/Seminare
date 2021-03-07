import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.interpolate import interp1d
import os

pgf_with_latex = {
    "pgf.texsystem": "pdflatex",
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": [],
    "font.sans-serif": [],
    "font.monospace": [],
    "axes.labelsize": 10,
    "font.size": 14,
    "legend.fontsize": 12,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "pgf.preamble": [
        r"\usepackage[utf8x]{inputenc}",
        r"\usepackage[T1]{fontenc}",
        r"\usepackage{amsmath}",
        ]
    }
mpl.rcParams.update(pgf_with_latex)

def plot_cantor_set(n = 4, show = False, save = False, path = './'):
    # Cantor-Menge, https://stackoverflow.com/questions/49071229/how-to-draw-a-zoomable-cantor-set-from-a-python-list

    line = [0,1]
    depth = n
    fig = plt.figure(figsize = (6, 4))
    ax = fig.add_subplot(111) 
    
    def divide(line, level = 0):
        ax.plot(line,[level, level], color = "blue", lw = 5, solid_capstyle = "butt")
        if level < depth:
            s = np.linspace(line[0], line[1], 4)
            divide(s[:2], level + 1)
            divide(s[2:], level + 1)

    divide(line)
    fig.gca().invert_yaxis()
    levels = np.linspace(0, depth, depth + 1)
    ax.set_yticks(levels)
    labels = [r'$n = {{{}}}$'.format(int(levels[i])) for i in range(depth + 1)]
    ax.set_yticklabels(labels)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.yaxis.set_tick_params(length = 0)
    
    if show:
        plt.show()
    if save:
        plt.savefig(path + 'cantor_set_{}.pdf'.format(n))
        plt.close()

def cantor(n):
    return [0.] + cant(0., 1., n) + [1.]

def cant(x, y, n):
    # rekursive Cantor-Funktion, https://stackoverflow.com/questions/17809817/cantor-ternary-set-in-python-or-c/17810389#17810389
    
    if n == 0:
        return []

    new_pts = [2.*x/3. + y/3., x/3. + 2.*y/3.]
    return cant(x, new_pts[0], n-1) + new_pts + cant(new_pts[1], y, n-1)

def plot_cantor_function(n = 4, show = False, save = False, path = './'):
    x = np.array(cantor(n))
    y = np.linspace(0., 1., 2**n + 1)
    ploty = np.concatenate([[y[0]], np.array([y[1:-1],y[1:-1]]).T.flatten(), [y[-1]]])
    
    fig = plt.figure(figsize = (6, 4))
    ax = fig.add_subplot(111) 
    ax.plot(x, ploty, color = 'blue')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    if show:
        plt.show()
    if save:
        plt.savefig(path + 'cantor_f_{}.pdf'.format(n))
        plt.close()

def ex_func_seq(x, n):
    if x < -1. / n:
        return 1. + x
    elif x >= -1. / n and x <= 1. / n:
        return -n / 2. * x**2 + 1 - 1. / (2. * n)
    elif x > 1. / n:
        return 1. - x

def ex_func_lim(x):
    if x <= 0.:
        return 1. + x
    else:
        return 1. - x

def ex_func_seq_deriv(x, n):
    if x < -1. / n:
        return 1.
    elif x >= -1. / n and x <= 1. / n:
        return -n * x
    elif x > 1. / n:
        return -1.

def ex_func_lim_deriv(x):
    if x < 0.:
        return 1.
    elif x == 0.:
        return 0
    elif x > 0.:
        return -1.

def plot_seq(n = 3, show = False, save = False, path = './'):
    nx = 500
    x = np.linspace(-1., 1., nx)
    un, u = np.zeros(nx), np.zeros(nx)
    for i in range(nx):
        un[i] = ex_func_seq(x[i], n)
        u[i] = ex_func_lim(x[i])
    
    fig = plt.figure(figsize = (6, 4))
    ax = fig.add_subplot(111) 
    ax.plot(x, u, color = 'red', label = r'$u$')
    ax.plot(x, un, color = 'blue', label = r'$u_n$')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xticks([-1., -1. / n, 0., 1. / n, 1.])
    ax.set_xticklabels([r'$-1$', r'$- 1 / n$', r'$0$', r'$+ 1 / n$', r'$1$'])
    ax.set_yticks([0., 1 - 1. / (2 * n), 1.])
    ax.set_yticklabels([r'$0$', r'$1 - \frac{1}{2 n}$', r'$1$'])
    plt.legend(loc = 'best')
    
    if show:
        plt.show()
    if save:
        plt.savefig(path + 'un.pdf')
        plt.close()

def plot_seq_deriv(n = 3, show = False, save = False, path = './'):
    nx = 500
    x = np.linspace(-1., 1., nx)
    unp, up = np.zeros(nx), np.zeros(nx)
    for i in range(nx):
        unp[i] = ex_func_seq_deriv(x[i], n)
        up[i] = ex_func_lim_deriv(x[i])
    
    fig = plt.figure(figsize = (6, 4))
    ax = fig.add_subplot(111) 
    ax.plot(x, up, color = 'red', label = r'$u^\prime$')
    ax.plot(x, unp, color = 'blue', label = r'$u_n^\prime$')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xticks([-1., -1. / n, 0., 1. / n, 1.])
    ax.set_xticklabels([r'$-1$', r'$- 1 / n$', r'$0$', r'$+ 1 / n$', r'$1$'])
    ax.set_yticks([-1., 0., +1.])
    ax.set_yticklabels([r'$-1$', r'$0$', r'$1$'])
    plt.legend(loc = 'best')
    
    if show:
        plt.show()
    if save:
        plt.savefig(path + 'unp.pdf')
        plt.close()

def plot_extension_u(show = False, save = False, path = './'):
    
    a = -1.
    b = 1.
    ap = -0.7
    bp = 0.4
    
    # interpolate between some arbitrary points for initial function u on [a,b]
    nx = 1000
    xab = np.linspace(a, b, nx)
    xp = [a, 0.2 * a + 0.8 * b, 0.5 * a + 0.5 * b, 0.25 * a + 0.75 * b, b]
    yp = [1., 0.2, 0.2, 0.5, 0.8]
    
    u1 = interp1d(xp[:-1], yp[:-1], kind = 'cubic')
    u2 = interp1d(xp[-2:], yp[-2:], kind = 'linear')
    
    def u(x):
        if x < a or x > b:
            return 0.
        elif x < xp[-2]:
            return u1(x)
        elif x >= xp[-2]:
            return u2(x)
            
    u_arr = np.zeros(nx)
    for i in range(nx):
        u_arr[i] = u(xab[i])
    
    
    # choose eta as cubic between ap and bp
    
    mtr = np.array([[ap**3, ap**2, ap, 1.],[bp**3, bp**2, bp, 1.], [3. * ap**2, 2. * ap, 1., 0.], [3. * bp**2, 2. * bp, 1., 0.]])
    alpha = np.linalg.inv(mtr) @ np.array([1., 0., 0., 0.])
    
    def eta(x):
        if x < ap:
            return 1.
        elif x > bp:
            return 0.
        else:
            return alpha[0] * x**3 + alpha[1] * x**2 + alpha[2] * x + alpha[3]
    
    nx_ext = 3 * nx
    x_extend = np.linspace(3. * a, 3. * b, nx_ext)
    eta_arr = np.zeros(nx_ext)
    for i in range(nx_ext):
        eta_arr[i] = eta(x_extend[i])
    
    def U1(x):
        if x >= a:
            return eta(x) * u(x)
        else:
            mirr = a + (a - x)
            return U1(mirr)
    
    def U2(x):
        if x <= b:
            return (1. - eta(x)) * u(x)
        else:
            mirr = b - (x - b)
            return U2(mirr)
    
    U1_arr, U2_arr = np.zeros(nx_ext), np.zeros(nx_ext)
    for i in range(nx_ext):
        U1_arr[i] = U1(x_extend[i])
        U2_arr[i] = U2(x_extend[i])
    
    U_arr = U1_arr + U2_arr
    
    def mollifier(x, eps):
        if abs(x) < eps:
            return np.exp(1. / ((x / eps)**2 - 1.))
        else: return 0
    
    eps = 0.2
    phi = np.zeros(nx_ext)
    for i in range(nx_ext):
        phi[i] = mollifier(x_extend[i], eps)
    
    phi = phi / (np.sum(phi) * (x_extend[1] - x_extend[0]))
    U_conv = np.fft.fftshift(np.fft.irfft(np.fft.rfft(U_arr) * np.fft.rfft(phi)) * (x_extend[1] - x_extend[0]))
    
    # initial function plot
    
    fig = plt.figure(figsize = (6, 4))
    ax = fig.add_subplot(111) 
    
    vwidth = 0.6
    ax.axvline(x = a, color = 'black', linestyle = '--', linewidth = vwidth)
    ax.axvline(x = b, color = 'black', linestyle = '--', linewidth = vwidth)
    ax.axvline(x = ap, color = 'black', linestyle = '--', linewidth = vwidth)
    ax.axvline(x = bp, color = 'black', linestyle = '--', linewidth = vwidth)
    
    ax.plot(xab, u_arr, color = 'blue', label = r'$u$')
    ax.set_xlim(1.2 * a, 1.2 * b)
    ax.set_ylim(-1,1.5)
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_yticks([0.])
    
    ax.set_xticks([a, ap, bp, b])
    ax.set_xticklabels([r'$a$', r'$a^\prime$', r'$b^\prime$', r'$b$'])
    plt.legend(loc = 'best')
    
    if show:
        plt.show()
    if save:
        plt.savefig(path + 'ext_init.pdf')
        plt.close()
    
    # eta plot
    
    fig = plt.figure(figsize = (6, 4))
    ax = fig.add_subplot(111) 
    
    vwidth = 0.6
    ax.axvline(x = a, color = 'black', linestyle = '--', linewidth = vwidth)
    ax.axvline(x = b, color = 'black', linestyle = '--', linewidth = vwidth)
    ax.axvline(x = ap, color = 'black', linestyle = '--', linewidth = vwidth)
    ax.axvline(x = bp, color = 'black', linestyle = '--', linewidth = vwidth)
    
    ax.plot(x_extend, eta_arr, color = 'green', label = r'$\eta$')
    ax.plot(x_extend, 1 - eta_arr, color = 'brown', label = r'$1 - \eta$')
    ax.set_xlim(1.5 * a, 1.5 * b)
    ax.set_ylim(0,1.2)
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_yticks([0., 1.])
    
    ax.set_xticks([a, ap, bp, b])
    ax.set_xticklabels([r'$a$', r'$a^\prime$', r'$b^\prime$', r'$b$'])
    plt.legend(loc = 'best')
    
    if show:
        plt.show()
    if save:
        plt.savefig(path + 'ext_eta.pdf')
        plt.close()
    
    # U_1, U_2 plot
    
    fig = plt.figure(figsize = (6, 4))
    ax = fig.add_subplot(111) 
    
    vwidth = 0.6
    ax.axvline(x = a, color = 'black', linestyle = '--', linewidth = vwidth)
    ax.axvline(x = b, color = 'black', linestyle = '--', linewidth = vwidth)
    ax.axvline(x = ap, color = 'black', linestyle = '--', linewidth = vwidth)
    ax.axvline(x = bp, color = 'black', linestyle = '--', linewidth = vwidth)
    
    ax.plot(x_extend, U1_arr, color = 'green', label = r'$U_1$')
    ax.plot(x_extend, U2_arr, color = 'brown', label = r'$U_2$')
    ax.set_xlim(3 * a, 3 * b)
    ax.set_ylim(-1,1.5)
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_yticks([0.])
    
    ax.set_xticks([a, ap, bp, b])
    ax.set_xticklabels([r'$a$', r'$a^\prime$', r'$b^\prime$', r'$b$'])
    plt.legend(loc = 'best')
    
    if show:
        plt.show()
    if save:
        plt.savefig(path + 'ext_mirror.pdf')
        plt.close()
    
    # finished
    
    fig = plt.figure(figsize = (6, 4))
    ax = fig.add_subplot(111)
    
    vwidth = 0.6
    ax.axvline(x = a, color = 'black', linestyle = '--', linewidth = vwidth)
    ax.axvline(x = b, color = 'black', linestyle = '--', linewidth = vwidth)
    ax.axvline(x = ap, color = 'black', linestyle = '--', linewidth = vwidth)
    ax.axvline(x = bp, color = 'black', linestyle = '--', linewidth = vwidth) 
    
    ax.plot(x_extend, U_arr, color = 'red', label = r'$U$')
    ax.plot(x_extend, U_conv, color = 'orange', label = r'${\cal S}_\varepsilon U$', alpha = 0.9)
    ax.plot(xab, u_arr, color = 'blue', label = r'$u$')
    ax.set_xlim(3 * a, 3 * b)
    ax.set_ylim(-1,1.5)
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_yticks([0.])
    
    ax.set_xticks([a, ap, bp, b])
    ax.set_xticklabels([r'$a$', r'$a^\prime$', r'$b^\prime$', r'$b$'])
    plt.legend(loc = 'best')
    
    if show:
        plt.show()
    if save:
        plt.savefig(path + 'ext_final.pdf')
        plt.close()
    

if __name__ == '__main__':
    
    folder = './figures/'
    
    if not os.path.isdir(folder):
        os.makedirs(folder)
    
    plot_cantor_set(save = True, path = folder)
    for i in range(3):
        plot_cantor_function(n = i, save = True, path = folder)
    
    plot_cantor_function(n = 10, save = True, path = folder)
    plot_seq(save = True, path = folder)
    plot_seq_deriv(save = True, path = folder)
    plot_extension_u(save = True, path = folder)
   
