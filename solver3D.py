from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
print('numpy: '+np.version.full_version)
import matplotlib.animation as animation
import matplotlib
print('matplotlib: '+matplotlib.__version__)



def SWAP(x0, x):
    return x, x0


def add_source(N, x, s, dt):
    for i in range(1, N):
        for j in range(1, N):
            for k in range(1, N):
                x[i][j][k] = dt* s[i][j][k]


def set_bnd(N, b, x):
    for i in range(1, N):
        for z in range(1, N):
            x[0][ i][ z] = -x[1][ i][ z] if b == 1 else x[1][ i][ z]
            x[N + 1][ i][ z] = -x[N][ i][ z] if b == 1 else x[N][ i][ z]
            x[i][ 0][ z] = -x[i][ 1][ z] if b == 2 else x[i][ 1][ z]
            x[i][ N + 1][ z] = -x[i][ N][ z] if b == 2 else x[i][ N][ z]
            x[i][ z][ 0] = -x[i][ z][ 1] if b == 3 else x[i][ z][ 1]
            x[i][ z][ N + 1] = -x[i][ z][ N] if b == 3 else x[i][ z][ N]

    x[0][ 0][ 0] = 1.0 / 3 * (x[1][ 0][ 0] + x[0][ 1][ 0] + x[0][ 0][ 1])
    x[N + 1][ 0][ 0] = 1.0 / 3 * (x[N][ 0][ 0] + x[N + 1][ 1][ 0] + x[N + 1][ 0][ 1])
    x[0][ N + 1][ 0] = 1.0 / 3 * (x[0][ N][ 0] + x[1][ N + 1][ 0] + x[0][ N + 1][ 1])
    x[0][ 0][ N + 1] = 1.0 / 3 * (x[0][ 0][ N] + x[1][ 0][ N + 1] + x[0][ 1][ N + 1])
    x[N + 1][ N + 1][ 0] = 1.0 / 3 * (x[N][ N + 1][ 0] + x[N + 1][ N][ 0] + x[N + 1][ N + 1][ 1])
    x[N + 1][ 0][ N + 1] = 1.0 / 3 * (x[N][ 0][ N + 1] + x[N + 1][ 1][ N + 1] + x[N + 1][ 0][ N])
    x[0][ N + 1][ N + 1] = 1.0 / 3 * (x[1][ N + 1][ N + 1] + x[0][ N][ N + 1] + x[0][ N + 1][ N])
    x[N + 1][ N + 1][ N + 1] = 1.0 / 3 * (x[N][ N + 1][ N + 1] + x[N + 1][ N][ N + 1] + x[N + 1][ N + 1][ N])


def lin_solve(N, b, x, x0, a, c):
    for m in range(0, 20):
        for i in range(1, N):
            for j in range(1, N):
                for k in range(1, N):
                    x[i][ j][ k] = (x0[i][ j][ k] + a * (
                            x[i - 1][ j][ k] + x[i + 1][ j][ k] + x[i][ j - 1][ k] + x[i][ j + 1][ k] + x[i][ j][ k + 1] + x[
                        i][ j][ k - 1])) / c
        set_bnd(N, b, x)


def diffuse(N, b, x, x0, diff, dt):
    a = dt * diff * N * N * N
    lin_solve(N, b, x, x0, a, 1 + 6 * a)


def advect(N, b, d, d0, u, v, w, dt):
    dt0 = dt * N
    for i in range(1, N):
        for j in range(1, N):
            for k in range(1, N):
                x = i - dt0 * u[i][ j][ k]
                y = j - dt0 * v[i][ j][ k]
                z = k - dt0 * w[i][ j][ k]

                if x < 0.5: x = 0.5
                if x > N + 0.5: x = N + 0.5
                i1 = int(x)
                i2 = i1 + 1
                if y < 0.5: y = 0.5
                if y > N + 0.5: y = N + 0.5
                j1 = int(y)
                j2 = j1 + 1
                if z < 0.5: z = 0.5
                if z > N + 0.5: z = N + 0.5
                k1 = int(y)
                k2 = k1 + 1

                d[i][ j][ k] = d0[i1][ j1][ k1] * (i2 - x) * (j2 - y) * (k2 - z) + d0[i2][ j1][ k1] * (x - i1) * (j2 - y) * (
                        k2 - z) \
                             + d0[i2][ j1][ k2] * (x - i1) * (j2 - y) * (z - k1) + d0[i2][ j2][ k1] * (x - i1) * (
                                     y - j1) * (k2 - z) \
                             + d0[i2][ j2][ k2] * (x - i1) * (y - j1) * (z - k1)
    set_bnd(N, b, d)


def project(N, u, v, w, p, div):
    for i in range(1, N):
        for j in range(1, N):
            for k in range(1, N):
                div[i][ j][ k] = -0.5 * (u[i + 1][ j][ k] - u[i - 1][ j][ k] + v[i][ j + 1][ k] - v[i][ j - 1][ k]
                                       + w[i][ j][ k + 1] - w[i][ j][ k + 1]) / N
                p[i][ j][ k] = 0
    set_bnd(N, 0, div)
    set_bnd(N, 0, p)

    lin_solve(N, 0, p, div, 1, 6)

    for i in range(1, N):
        for j in range(1, N):
            for k in range(1, N):
                u[i][ j][ k] -= 0.5 * N * (p[i + 1][ j][ k] - p[i - 1][ j][ k])
                v[i][ j][ k] -= 0.5 * N * (p[i][ j + 1][ k] - p[i][ j - 1][ k])
                w[i][ j][ k] -= 0.5 * N * (p[i][ j][ k + 1] - p[i][ j][ k - 1])

    set_bnd(N, 1, u)
    set_bnd(N, 2, v)


def dens_step(N, x, x0, u, v, w, diff, dt):
    add_source(N, x, x0, dt)
    x0, x = SWAP(x0, x)
    diffuse(N, 0, x, x0, diff, dt)
    x0, x = SWAP(x0, x)
    advect(N, 0, x, x0, u, v, w, dt)


def vel_step(N, u, v, w, u0, v0, w0, visc, dt):
    add_source(N, u, u0, dt)
    add_source(N, v, v0, dt)
    add_source(N, w, w0, dt)
    u0, u = SWAP(u0, u)
    diffuse(N, 1, u, u0, visc, dt)
    v0, v = SWAP(v0, v)
    diffuse(N, 2, v, v0, visc, dt)
    w0, w = SWAP(w0, w)
    diffuse(N, 2, w, w0, visc, dt)
    project(N, u, v, w, u0, v0)
    u0, u = SWAP(u0, u)
    v0, v = SWAP(v0, v)
    w0, w = SWAP(w0, w)
    advect(N, 1, u, u0, u0, v0, w0, dt)
    advect(N, 2, v, v0, u0, v0, w0, dt)
    advect(N, 3, w, v0, u0, v0, w0, dt)
    project(N, u, v, w, u0, v0)





#
#
#
# Nfrm = 10
# fps = 10
#
# def animation():
#
#     N = 64
#     dt = 0.1
#     diff = 0.0
#     visc = 0.0
#     force = 5.0
#     source = 100.0
#     print( "Using defaults : N= " , N, " dt= ", dt, " diff= ", diff, " visc= ", visc, " force = ", force, " source= ", source)
#     # Set up Figure and 3D Axes
#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')
#     # Set axes limits so that the whole image is included
#     ax.set(xlim=(0, 70), ylim=(0, 70), zlim=(0, 70))
#     # Draw a blank line
#     line, = ax.plot([], [])
#
#
#
#     # Define data
#     # x = np.linspace(0, 64, 1)
#     # y = np.linspace(0, 64, 1)
#     # z = np.linspace(0, 64, 1)
#     # ax.plot(x, y,z)
#
#     plt.show()
#
# if __name__ == '__main__':
#     animation()