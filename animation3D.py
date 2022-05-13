import matplotlib.pyplot as plt
import solver3D
import numpy as np
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D

def add_forces(u, v, w, xyz_forces):
    if xyz_forces[0] == -1: return
    xyz1 = [0,0,0]
    xyz1[0] = N if (xyz_forces[0]>N+1) else xyz_forces[0]+1
    xyz1[1] = N if (xyz_forces[0] > N + 1) else xyz_forces[1] + 1
    xyz1[2] = N if (xyz_forces[0] > N + 1) else xyz_forces[2] + 1
    u[int(xyz1[0])][ int(xyz1[1])][int(xyz1[2])] = force*(xyz1[0]-xyz_forces[0])
    v[int(xyz1[0])][int(xyz1[1])][ int(xyz1[2])] = force * (xyz1[1] - xyz_forces[1])
    w[int(xyz1[0])][int(xyz1[1])][ int(xyz1[2])] = force * (xyz1[2] - xyz_forces[2])


def add_density(d, xyz, source):
    if xyz[0] == -1: return
    for i in range(len(xyz)):
        d[xyz[0]][ xyz[1]][ xyz[2]] = source


def draw_density(N, h, dens):
    for i in range(1, N):
        x = (i - 0.5) * h
        for j in range(1, N):
            y = (j - 0.5) * h
            for k in range(1, N):
                z = (k - 0.5) * h
                d000 = dens[i][j][z]
                d001 = dens[i][j][z + 1]
                d010 = dens[i][ j+1][z]
                d011 = dens[i+1][j][z]
                d100 = dens[i][ j][z-1]
                d101 = dens[i][ j-1][z]
                d110 = dens[i-1][ j][z]
                plt.plot(x,y,z, color = d000)
                plt.plot(x, y, z+1, color=d001)
                plt.plot(x, y+1, z, color=d010)
                plt.plot(x+1, y, z, color=d011)
                plt.plot(x, y, z-1, color=d100)
                plt.plot(x, y-1, z, color=d101)
                plt.plot(x-1, y, z, color=d110)



def draw_velocity(N, h, u, v, w):

    for i in range(1, N):
        x = (i - 0.5) * h
        for j in range(1, N):
            y = (j - 0.5) * h
            for k in range(1, N):
                z = (k - 0.5) * h

                x_vector = [x, x+u[i][j][k]]
                y_vector = [y, y+v[i][j][k]]
                z_vector = [z, z+w[i][j][k]]
                ax.plot( x_vector, y_vector, z_vector, color = "blue")

def step_time(N, u, v, w, u_prev, v_prev, w_prev, dens, dens_prev, visc, diff, dt, xyz_sources, xyz_forces):
    add_forces( u_prev, v_prev, w_prev, xyz_forces)
    print("ok1")
    add_density(dens_prev, xyz_sources, source)
    print("ok2")
    solver3D.vel_step(N, u, v, w, u_prev, v_prev, w_prev, visc, dt)
    print("ok3")
    solver3D.dens_step(N, dens, dens_prev, u, v, w, diff, dt)
    print("ok4")



fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set(xlim=(0, 70), ylim=(0, 70), zlim=(0, 70))
line1, = ax.plot([], [], [], lw=3)

N = 30
h = 1
dt = 0.1
diff = 0.0
visc = 0.0
force = 5.0
source = 100.0
size = (N+2)**3
number_frames = 20
u = []
v = []
w = []
w_prev = []
v_prev = []
u_prev = []
dens = []
dens_prev = []
xyz_sources = np.random.randint(1,N, size=(number_frames,3))
xyz_forces = np.random.randint(1,N, size=(number_frames,3))
for i in range(number_frames):
    if i%3 != 0 :
        xyz_sources[i] = [-1,-1,-1]
        xyz_forces[i] = [-1,-1,-1]

w = [[[0 for i in range(N+2)]for j in range(N+2)] for k in range(N+2)]
u = [[[0 for i in range(N+2)]for j in range(N+2)] for k in range(N+2)]
v = [[[0 for i in range(N+2)]for j in range(N+2)] for k in range(N+2)]
w_prev = [[[0 for i in range(N+2)]for j in range(N+2)] for k in range(N+2)]
v_prev = [[[0 for i in range(N+2)]for j in range(N+2)] for k in range(N+2)]
u_prev = [[[0 for i in range(N+2)]for j in range(N+2)] for k in range(N+2)]
dens = [[[0 for i in range(N+2)]for j in range(N+2)] for k in range(N+2)]
dens_prev = [[[0 for i in range(N+2)]for j in range(N+2)] for k in range(N+2)]


def animate(i):

    step_time(N, u, v, w, u_prev, v_prev, w_prev, dens, dens_prev, visc, diff, dt, xyz_sources[i], xyz_forces[i])
    draw_velocity(N, h, u, v, w)
# anim1 = FuncAnimation(fig, animate,
#                       frames=number_frames, interval=50, blit=True)

animate(0)

plt.show()
