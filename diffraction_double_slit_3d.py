# Double slit diffraction in 3D
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import cm
import mpl_toolkits.mplot3d.art3d as art3d


def draw_double_slit(s_width, s_right_center, s_left_center):
    sw = s_width / 2
    x_wall_right = [0, 0, 0, 0]
    y_wall_right = [s_right_center + sw, y_max, y_max, s_right_center + sw]
    z_wall_right = [z_max, z_max, z_min, z_min]
    poly = list(zip(x_wall_right, y_wall_right, z_wall_right))
    ax.add_collection3d(art3d.Poly3DCollection([poly], color='gray'))

    x_wall_center = [0, 0, 0, 0]
    y_wall_center = [s_left_center + sw, s_right_center - sw, s_right_center - sw, s_left_center + sw]
    z_wall_center = [z_max, z_max, z_min, z_min]
    poly = list(zip(x_wall_center, y_wall_center, z_wall_center))
    ax.add_collection3d(art3d.Poly3DCollection([poly], color='gray'))

    x_wall_left = [0, 0, 0, 0]
    y_wall_left = [s_left_center - sw, y_min, y_min, s_left_center - sw]
    z_wall_left = [z_max, z_max, z_min, z_min]
    poly = list(zip(x_wall_left, y_wall_left, z_wall_left))
    ax.add_collection3d(art3d.Poly3DCollection([poly], color='gray'))


def set_axis():
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_zlim(z_min, z_max)
    ax.set_title('Double slit diffraction ; sin(k*x - omega*t) (step as t)')
    ax.set_xlabel('x * pi')
    ax.set_ylabel('y * pi')
    ax.set_zlabel('z')
    ax.set_box_aspect((1, 1, 0.5))
    ax.grid()


def update(f):
    ax.cla()  # Clear ax
    set_axis()

    global zz_superposed  # Global variables that changes in this function
    # Select mode (combination of parameters)
    mode = int(f / 20) % 6
    if mode == 0:
        k, omega, slit_width, slit_right_center, slit_left_center = parameters0
    elif mode == 1:
        k, omega, slit_width, slit_right_center, slit_left_center = parameters1
    elif mode == 2:
        k, omega, slit_width, slit_right_center, slit_left_center = parameters2
    elif mode == 3:
        k, omega, slit_width, slit_right_center, slit_left_center = parameters3
    elif mode == 4:
        k, omega, slit_width, slit_right_center, slit_left_center = parameters4
    elif mode == 5:
        k, omega, slit_width, slit_right_center, slit_left_center = parameters5
    else:
        k, omega, slit_width, slit_right_center, slit_left_center = parameters0

    # Draw text
    ax.text(x_min, y_min, z_max * 3.0, "Step=" + str(f) + ", Parameter-combination No." + str(mode))
    ax.text(x_min, y_min, z_max * 2.7, "k=" + str(k) + ", omega=" + str(omega))
    ax.text(x_min, y_min, z_max * 2.4, "Distance of double-slit=" + str(slit_right_center - slit_left_center))
    ax.text(x_min, y_min, z_max * 2.1, "Width of slits=" + str(slit_width))
    # Draw double slit
    draw_double_slit(slit_width, slit_right_center, slit_left_center)
    # Draw a wave, Note: np.pi in sine is for adjustment x axis as x * pi
    # Superpose(sum) sine value at xx,yy with distance from each super_position_points to xx,yy
    super_position_points = np.linspace(-slit_width / 2, slit_width / 2, num_of_superpose)
    zz_superposed = xx * 0. + yy * 0.
    for i in super_position_points:
        length = np.sqrt((xx - 0.) ** 2 + (yy - (i + slit_right_center)) ** 2)
        zz_superposed = zz_superposed + np.sin(k * length * np.pi - omega * f)
    for i in super_position_points:
        length = np.sqrt((xx - 0.) ** 2 + (yy - (i + slit_left_center)) ** 2)
        zz_superposed = zz_superposed + np.sin(k * length * np.pi - omega * f)
    zz_superposed = zz_superposed / len(super_position_points) / 2
    ax.plot_surface(xx, yy, zz_superposed, rstride=1, cstride=1, cmap=cm.coolwarm, alpha=0.6)


# Global variables
x_min = -0.5
x_max = 4.
y_min = -2.
y_max = 2.
z_min = -4.
z_max = 4.

# Prepare mesh grid
x = np.arange(0, x_max, 0.05)
y = np.arange(y_min, y_max, 0.05)
xx, yy = np.meshgrid(x, y)
zz_superposed = xx * 0. + yy * 0.

# Double slit initial parameters
slit_width0 = 0.4
slit_right_center0 = 0.4
slit_left_center0 = -0.4

# Wave initial parameters
k0 = 1.     # Wave number
omega0 = 1.  # Angular velocity

# Combination of parameters
parameters0 = [k0, omega0, slit_width0, slit_right_center0, slit_left_center0]
parameters1 = [k0 * 4, omega0, slit_width0, slit_right_center0, slit_left_center0]
parameters2 = [k0 * 8, omega0, slit_width0, slit_right_center0, slit_left_center0]
parameters3 = [k0, omega0, slit_width0 / 2., slit_right_center0, slit_left_center0]
parameters4 = [k0 * 4, omega0, slit_width0 / 2., slit_right_center0, slit_left_center0]
parameters5 = [k0 * 8, omega0, slit_width0 / 2., slit_right_center0, slit_left_center0]

# Number of superposition points in each slits
num_of_superpose = 9

# Generate figure and axes
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Draw animation
set_axis()
anim = animation.FuncAnimation(fig, update, interval=100)
plt.show()
