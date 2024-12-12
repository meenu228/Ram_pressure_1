# We calculate the elliptical orbits using leap frog method

import matplotlib

#matplotlib.use('TkAgg')

import matplotlib.pyplot as plt
from matplotlib import animation

from matplotlib.lines import Line2D

from numpy import sin, cos
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation

import astropy.units as u
from astropy.constants import G


matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
# plt.style.use('seaborn-talk')

matplotlib.rc('text', usetex=True)
matplotlib.rc('text', usetex=True)
matplotlib.rc('axes', linewidth=2)
matplotlib.rc('font', weight='bold')

plt.rc('font', family='serif')
plt.rc('xtick', labelsize='x-small')
plt.rc('ytick', labelsize='x-small')


plt.style.use('classic')
#matplotlib.rcParams['mathtext.rm'] = 'roman'

classic_changed = {
    'lines.linewidth': 1.0,
    'lines.dashed_pattern' : [6, 6],
    'lines.dashdot_pattern' : [4.8, 1.2, 0.8, 1.2],
    'lines.dotted_pattern' : [1, 3],
    'lines.scale_dashes': False}



G = G.to('Mpc3/(Msun Myr2)')

f_g = 0.8
O_b = 0.015 / (0.7 ** 2)
O_m = 0.3

omega_m = 0.316
omega_l = 0.684
omega_k = 0.0
w0 = -1.0
wa = 0.0
H0 = 0.00010227121650537077
f_uni = 0.158
h = 0.67

M_gal = 1.0e11 * u.Msun
M_ICM = 2.0e15 * u.Msun

z = 0.0

omega_m_z = (omega_m * ((1.0 + z) ** 3)) / (omega_m * ((1.0 + z) ** 3) + omega_k * ((1.0 + z) ** 2) + omega_l)
d = omega_m_z - 1.0
delta_c = 18.0 * np.pi ** 2 + 82.0 * d - 39.0 * d ** 2


def f_de(a):
    epsilon = 0.000000001
    return -3.0 * (1.0 + w0) + 3.0 * wa * ((a - 1.0) / np.log(a - epsilon) - 1.0)


def Esqr(a):
    return omega_m * pow(a, -3) + omega_k * pow(a, -2) + omega_l * pow(a, f_de(a))


def H(a):
    return h * H0 * np.sqrt(Esqr(a))


def f1(x, R_CORE, c1):
    beta = 0.6
    FUNCT1 = (x ** 2) / ((1 + x ** 2) ** (3.0 * beta / 2.0))
    return FUNCT1


def f2(x, R_vir, c1):
    c = c1
    FUNCT2 = ((x + 0.000001) ** 2) / ((c * (x + 0.000001)) * ((1.0 + c * (x + 0.000001)) ** 2))
    return FUNCT2


def M_NFW(r_c, R_200):
    H_z = H(1.0 / (1 + z)) * (1 / u.Myr)
    rho_c = 3.0 * H_z ** 2 / (8.0 * np.pi * (G.to('Mpc3/(Msun Myr2)')))
    c = 10
    delta_c = (200.0 / 3.0) * (c ** 3 / (np.log(c + 1.0) - (c / (1 + c))))
    ymax = r_c / R_200
    R = R_200.value
    M = 4 * np.pi * delta_c * rho_c * (R_200 ** 3) * integrate.romberg(f2, 0.0, ymax.value, args=(R, c), divmax=100)
    
    # print integrate.romberg(f2, 0.0, ymax.value,args=(R,c), divmax=100),rho_c*delta_c

    return M


def Beta_M(r_c, R_200):
    H_z = H(1.0 / (1 + z)) *  (1 / u.Myr)
    c = 10.0

    f_uni = 0.158
    rho_c = (3.0 * (H_z ** 2)) / (8.0 * np.pi * (G.to('Mpc3/(Msun Myr2)')))
    delta_c = (200.0 / 3.0) * ((c ** 3) / (np.log(c + 1.0) - (c / (1 + c))))

    R = R_200.value / 20.0
    R1 = R_200.value
    del1 = 1.0
    f_g = 0.8
    O_b = 0.022 / (0.7 ** 2)
    O_m = 0.316

    DELTA = f_g * O_b / O_m    # f_uni/(1.0 - f_uni)#del1*f_uni

    M_BG = integrate.romberg(f2, 0.0, 1.0, args=(R1, c), divmax=100)
    CONST = integrate.romberg(f1, 0.0, 20.0, args=(R, c), divmax=100)
    rho_0 = 0.5 * DELTA * (20.0 ** 3) * (M_BG / CONST) * (rho_c * delta_c)

    ymax = r_c / (R * u.Mpc)
    c1 = 10.0

    # print rho_c,delta_c

    M = 4.0 * np.pi * ((R_200 / 20.0) ** 3) * (rho_0) * integrate.romberg(f1, 0.0, ymax.value, args=(R, c1), divmax=100)

    return M


def dist(a, b, c, d):
    #print("meenu is",np.sqrt((a - b) ** 2 + (c - d) ** 2))
    return np.sqrt((a - b) ** 2 + (c - d) ** 2)


def animate(i):
    theta = np.pi / 6.0
    # setting up the initial conditions
    r_0 = ((1.0 - e) / (1.0 + q)) * a.value
    v_0 = ((1.0 / (1.0 + q)) * np.sqrt((1.0 + e) / (1.0 - e)) * np.sqrt(G * (M_ICM + M_gal) / a)).to('Mpc/Myr')
    r_0 = 2 * R_vir.value
    v_0 = 0.5 * V200
    #print("V200 is",V200)

    if (i == 0):
        x1[i, 0] = r_0 * np.cos(theta)  # -400.011995319#
        y1[i, 0] = r_0 * np.sin(theta)  # 363.865867513#

        x2[i, 0] = 0.0  # 0.627784919989#
        y2[i, 0] = 0.0  # -0.571056636293#

        vx1[i, 0] = -v_0.value  # 86.5591029587#
        vy1[i, 0] = 0.0  # v_0.value  # -28.7996765414#

        vx2[i, 0] = 0.0  # -0.135847174988#
        vy2[i, 0] = 0.0  # 0.0451986511525#

        ax1[i, 0] = (G.value * M_ICM.value / (dist(x1[i, 0], x2[i, 0], y1[i, 0], y2[i, 0]) ** 2)) * (
        (x2[i, 0] - x1[i, 0]) / (dist(x1[i, 0], x2[i, 0], y1[i, 0], y2[i, 0])))

        ay1[i, 0] = (G.value * M_ICM.value / (dist(x1[i, 0], x2[i, 0], y1[i, 0], y2[i, 0]) ** 2)) * (
        (y2[i, 0] - y1[i, 0]) / (dist(x1[i, 0], x2[i, 0], y1[i, 0], y2[i, 0])))

        ax2[i, 0] = (G.value * M_gal.value / (dist(x1[i, 0], x2[i, 0], y1[i, 0], y2[i, 0]) ** 2)) * (
        (x1[i, 0] - x2[i, 0]) / (dist(x1[i, 0], x2[i, 0], y1[i, 0], y2[i, 0])))

        ay2[i, 0] = (G.value * M_gal.value / (dist(x1[i, 0], x2[i, 0], y1[i, 0], y2[i, 0]) ** 2)) * (
        (y1[i, 0] - y2[i, 0]) / (dist(x1[i, 0], x2[i, 0], y1[i, 0], y2[i, 0])))

        xdata.append(x1[i, 0])
        ydata.append(y1[i, 0])

    if (0 < i <= nstep - 1):
        # print x1[i-1,0],y1[i-1,0],x2[i-1,0],y2[i-1,0],vx1[i-1,0],vy1[i-1,0],vx2[i-1,0],vy2[i-1,0],ax1[i-1,0],ay1[i-1,0],ax2[i-1,0],ay2[i-1,0]


        vx1[i, 0] = vx1[i - 1, 0] + 0.5 * ax1[i - 1, 0] * dt
        vx2[i, 0] = vx2[i - 1, 0] + 0.5 * ax2[i - 1, 0] * dt
        vy1[i, 0] = vy1[i - 1, 0] + 0.5 * ay1[i - 1, 0] * dt
        vy2[i, 0] = vy2[i - 1, 0] + 0.5 * ay2[i - 1, 0] * dt

        x1[i, 0] = x1[i - 1, 0] + dt * vx1[i, 0]
        y1[i, 0] = y1[i - 1, 0] + dt * vy1[i, 0]
        x2[i, 0] = x2[i - 1, 0] + dt * vx2[i, 0]
        y2[i, 0] = y2[i - 1, 0] + dt * vy2[i, 0]

        M_r = M_NFW(dist(x1[i, 0], x2[i, 0], y1[i, 0], y2[i, 0]), R_200) + Beta_M(
            dist(x1[i, 0], x2[i, 0], y1[i, 0], y2[i, 0]), R_200)

        ax1[i, 0] = (G.value * M_r.value / (dist(x1[i, 0], x2[i, 0], y1[i, 0], y2[i, 0]) ** 2)) * (
        (x2[i, 0] - x1[i, 0]) / (dist(x1[i, 0], x2[i, 0], y1[i, 0], y2[i, 0])))

        ay1[i, 0] = (G.value * M_r.value / (dist(x1[i, 0], x2[i, 0], y1[i, 0], y2[i, 0]) ** 2)) * (
        (y2[i, 0] - y1[i, 0]) / (dist(x1[i, 0], x2[i, 0], y1[i, 0], y2[i, 0])))

        ax2[i, 0] = (G.value * M_gal.value / (dist(x1[i, 0], x2[i, 0], y1[i, 0], y2[i, 0]) ** 2)) * (
        (x1[i, 0] - x2[i, 0]) / (dist(x1[i, 0], x2[i, 0], y1[i, 0], y2[i, 0])))

        ay2[i, 0] = (G.value * M_gal.value / (dist(x1[i, 0], x2[i, 0], y1[i, 0], y2[i, 0]) ** 2)) * (
        (y1[i, 0] - y2[i, 0]) / (dist(x1[i, 0], x2[i, 0], y1[i, 0], y2[i, 0])))

        vx1[i, 0] = vx1[i, 0] + 0.5 * ax1[i, 0] * dt
        vx2[i, 0] = vx2[i, 0] + 0.5 * ax2[i, 0] * dt
        vy1[i, 0] = vy1[i, 0] + 0.5 * ay1[i, 0] * dt
        vy2[i, 0] = vy2[i, 0] + 0.5 * ay2[i, 0] * dt

        xdata.append(x1[i, 0])
        ydata.append(y1[i, 0])

        #print (x1[i, 0],y1[i, 0])

        '''

        if (i % 500.0 == 0):
            ax.arrow(x1[i, 0], y1[i, 0], 100 * vx1[i, 0] / (v_0.value), 0.0, head_width=10.0, head_length=10.0, fc='k',
                     ec='k')
            ax.arrow(x1[i, 0], y1[i, 0], 0.0, 100 * vy1[i, 0] / (v_0.value), head_width=10.0, head_length=10.0, fc='k',
                     ec='k')
        '''

    # line.set_data(xdata, ydata)


    scat.set_data(xdata, ydata)
    time_text.set_text('{:7s} {:3s} {:2s}'.format('time = ', str(i * dt / 1000)[:4], ' Gyr'))

    return scat, time_text


R_200 = (0.784 * ((M_ICM * h / (1.0e+8 * u.Msun)) ** (1.0 / 3.0)) * (
(omega_m_z * 18.0 * np.pi ** 2 / (omega_m * delta_c)) ** (1.0 / 3.0)) * (10.0 / (h * (1.0 + z))) * u.kpc).to('Mpc')

V200 = (23.4 * ((M_ICM * h / (1.0e8 * u.Msun)) ** (1.0 / 3.0)) * (
((omega_m * delta_c) / (omega_m_z * 18.0 * (np.pi ** 2))) ** (1.0 / 6.0)) * (((1.0 + z) / 10) ** (1.0 / 2.0)) * (
        u.km / u.second)).to('Mpc/Myr')

fig, ax = plt.subplots(figsize=(8, 8))
ax.grid()
time_text = ax.text(0.95, 0.95, '', horizontalalignment='right', verticalalignment='top', transform=ax.transAxes)
xdata, ydata = [], []
# line, = ax.plot([], [], lw=2.0,label=r'$\theta = \pi/4, V_{in} = 1 \times V_{200}$')

scat, = ax.plot([], [], linestyle='', ms=2, marker='o', color='b', label=r'$\mathrm{\theta = \pi/6, V_{in} = 0.5 \times V_{vir}}$')


def init():
    # line.set_data([], [])
    scat, = ax.plot([], [], linestyle='', ms=2, marker='o', color='b')
    ax.set_xlabel('x(Mpc)',  fontsize=16)
    ax.set_ylabel('y(Mpc)',  fontsize=16)
    ax.minorticks_on()
    ax.tick_params(axis='both', which='minor', length=5, width=2, labelsize=14)
    ax.tick_params(axis='both', which='major', length=8, width=2, labelsize=14)

    ax.set_xlim(-int(4 * R_200.value), int(4 * R_200.value))
    ax.set_ylim(-int(4 * R_200.value), int(4 * R_200.value))
    circle = plt.Circle((0.0, 0.0), 2 * R_vir.value, color='k', alpha=0.1)
    ax.add_artist(circle)
    ax.legend(loc=2)
    plt.annotate(s='', xy=(2 * R_200.value, 0), xytext=(0, 0), arrowprops=dict(arrowstyle='<->'))
    ax.text(R_vir.value / 2, -1.0, r'$\mathrm{R=2R_{vir}}$', fontsize=10)
    time_text.set_text('')

    return scat,


R_vir = R_200

q = (M_gal / M_ICM).value

e = 0.5

a = R_vir

tau = np.sqrt((4.0 * (np.pi ** 2) * a ** 3) / (G * (M_gal + M_ICM))).to('Myr').value
nstep = 5000
dt = tau / nstep

# dt=0.000000000001

x1 = np.zeros((nstep, 1))
print("x1 is", x1)
y1 = np.zeros((nstep, 1))
vx1 = np.zeros((nstep, 1))
vy1 = np.zeros((nstep, 1))
ax1 = np.zeros((nstep, 1))
ay1 = np.zeros((nstep, 1))

x2 = np.zeros((nstep, 1))
y2 = np.zeros((nstep, 1))
vx2 = np.zeros((nstep, 1))
vy2 = np.zeros((nstep, 1))
ay2 = np.zeros((nstep, 1))
ax2 = np.zeros((nstep, 1))

time_value = np.linspace(0, nstep, nstep)

X = time_value.astype(int)

ax.relim()  # reset intern limits of the current axes
ax.autoscale_view(True, True, True)

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=X, interval=20, blit=True)

#Writer = animation.writers['ffmpeg']
#writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

#anim.save('testing.gif', dpi=80, writer='imagemagick')

# anim.save('animation.mp4', writer=writer)

plt.show()

