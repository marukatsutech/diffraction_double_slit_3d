# Double slit diffraction in 3D (without clear ax)(embedded in tkinter)
import numpy as np
from matplotlib.figure import Figure
import matplotlib.animation as animation
from matplotlib import cm
import mpl_toolkits.mplot3d.art3d as art3d
import tkinter as tk
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)


def update_double_slit(s_width, s_center_distance):
    global wall_right, wall_center, wall_left
    s_right_center = s_center_distance / 2.
    s_left_center = - s_center_distance / 2.
    sw = s_width / 2.
    x_wall_right = [0, 0, 0, 0]
    y_wall_right = [s_right_center + sw, y_max, y_max, s_right_center + sw]
    z_wall_right = [z_max, z_max, z_min, z_min]
    poly = list(zip(x_wall_right, y_wall_right, z_wall_right))
    wall_right.remove()
    wall_right = ax.add_collection3d(art3d.Poly3DCollection([poly], color='gray'))

    x_wall_center = [0, 0, 0, 0]
    y_wall_center = [s_left_center + sw, s_right_center - sw, s_right_center - sw, s_left_center + sw]
    z_wall_center = [z_max, z_max, z_min, z_min]
    poly = list(zip(x_wall_center, y_wall_center, z_wall_center))
    wall_center.remove()
    wall_center = ax.add_collection3d(art3d.Poly3DCollection([poly], color='gray'))

    x_wall_left = [0, 0, 0, 0]
    y_wall_left = [s_left_center - sw, y_min, y_min, s_left_center - sw]
    z_wall_left = [z_max, z_max, z_min, z_min]
    poly = list(zip(x_wall_left, y_wall_left, z_wall_left))
    wall_left.remove()
    wall_left = ax.add_collection3d(art3d.Poly3DCollection([poly], color='gray'))


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


def reset():
    global is_play, cnt, surface1, zz_superposed
    is_play = False
    cnt = 0
    txt_step.set_text("Step=" + str(cnt))
    surface1.remove()
    zz_superposed = xx * 0. + yy * 0.
    surface1 = ax.plot_surface(xx, yy, zz_superposed, rstride=1, cstride=1, cmap=cm.coolwarm, alpha=0.6)


def change_slit_distance(value):
    global slit_center_distance
    reset()
    slit_center_distance_check = float(value)
    if slit_center_distance_check > slit_width:
        slit_center_distance = float(value)
        update_double_slit(slit_width, slit_center_distance)


def change_slit_width(value):
    global slit_width
    reset()
    slit_width_check = float(value)
    if slit_center_distance > slit_width_check:
        slit_width = float(value)
        update_double_slit(slit_width, slit_center_distance)


def change_omega(om):
    global omega
    reset()
    omega = float(om)
    txt_k_omega.set_text("k=" + str(k) + ", omega=" + str(omega))


def change_k(value):
    global k
    reset()
    k = float(value)
    txt_k_omega.set_text("k=" + str(k) + ", omega=" + str(omega))


def switch():
    global is_play
    if is_play:
        is_play = False
    else:
        is_play = True


def update(f):
    global cnt, zz_superposed, txt_step, txt_k_omega, txt_dist_slit, txt_width_slit, surface1
    slit_right_center = slit_center_distance / 2.
    slit_left_center = - slit_center_distance / 2.
    if is_play:
        # update text
        txt_step.set_text("Step=" + str(cnt))
        txt_k_omega.set_text("k=" + str(k) + ", omega=" + str(omega))
        txt_dist_slit.set_text("Distance of double-slit=" + str(slit_center_distance))
        txt_width_slit.set_text("Width of slits=" + str(slit_width))
        # Draw double slit
        update_double_slit(slit_width, slit_center_distance)
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
        surface1.remove()
        surface1 = ax.plot_surface(xx, yy, zz_superposed, rstride=1, cstride=1, cmap=cm.coolwarm, alpha=0.6)
        cnt += 1


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
slit_width_init = 0.4
slit_width = slit_width_init
slit_center_distance_init = 0.8
slit_center_distance = slit_center_distance_init

# Wave initial parameters
k_init = 1.     # Wave number
k = k_init
omega_init = 0.1  # Angular velocity
omega = omega_init

# For UI
is_play = False
cnt = 0

# Number of superposition points in each slits
num_of_superpose = 9

# Generate figure and axes
fig = Figure()
ax = fig.add_subplot(111, projection='3d')
set_axis()

# Generate items
txt_step = ax.text(x_min, y_min, z_max * 3.0, "Step=" + str(0))
txt_k_omega = ax.text(x_min, y_min, z_max * 2.7, "k=" + str(k_init) + ", omega=" + str(omega_init))
txt_dist_slit = ax.text(x_min, y_min, z_max * 2.4, "Distance of double-slit="
                        + str(slit_center_distance))
txt_width_slit = ax.text(x_min, y_min, z_max * 2.1, "Width of slits=" + str(slit_width_init))
surface1 = ax.plot_surface(xx, yy, zz_superposed, rstride=1, cstride=1, cmap=cm.coolwarm, alpha=0.6)

poly_right = list(zip([0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]))
wall_right = ax.add_collection3d(art3d.Poly3DCollection([poly_right], color='gray'))
poly_center = list(zip([0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]))
wall_center = ax.add_collection3d(art3d.Poly3DCollection([poly_center], color='gray'))
poly_left = list(zip([0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]))
wall_left = ax.add_collection3d(art3d.Poly3DCollection([poly_left], color='gray'))
update_double_slit(slit_width_init, slit_center_distance)

# Tkinter
root = tk.Tk()
root.title("Double slit diffraction")
canvas = FigureCanvasTkAgg(fig, root)
canvas.get_tk_widget().pack(expand=True, fill='both')

toolbar = NavigationToolbar2Tk(canvas, root)
canvas.get_tk_widget().pack()

# Play and pause button
btn_pp = tk.Button(root, text="Play/Pause", command=switch)
btn_pp.pack(side='left')

# Clear button
btn_clr = tk.Button(root, text="Reset", command=reset)
btn_clr.pack(side='left')

# Label and spinbox for k
lbl_k = tk.Label(root, text="k")
lbl_k.pack(side='left')
var_k = tk.StringVar(root)  # variable for spinbox-value
var_k.set(k_init)  # Initial value
spn_k = tk.Spinbox(
    root, textvariable=var_k, format="%.1f", from_=1., to=10., increment=1.,
    command=lambda: change_k(var_k.get()), width=5
    )
spn_k.pack(side='left')

# Label and spinbox for omega
lbl_omega = tk.Label(root, text=", Omega")
lbl_omega.pack(side='left')
var_omega = tk.StringVar(root)  # variable for spinbox-value
var_omega.set(k_init)  # Initial value
spn_omega = tk.Spinbox(
    root, textvariable=var_omega, format="%.1f", from_=0.1, to=2., increment=0.1,
    command=lambda: change_omega(var_omega.get()), width=5
    )
spn_omega.pack(side='left')

# Label and spinbox for the width of slits
lbl_ws = tk.Label(root, text=" ,Width of slits")
lbl_ws.pack(side='left')
var_ws = tk.StringVar(root)  # variable for spinbox-value
var_ws.set(slit_width_init)  # Initial value
spn_ws = tk.Spinbox(
    root, textvariable=var_ws, format="%.2f", from_=0.1, to=1.0, increment=0.1,
    command=lambda: change_slit_width(var_ws.get()), width=5
    )
spn_ws.pack(side='left')

# Label and spinbox for the distance of slits
lbl_ds = tk.Label(root, text=" ,Distance of slits")
lbl_ds.pack(side='left')
var_ds = tk.StringVar(root)  # variable for spinbox-value
var_ds.set(slit_center_distance_init)  # Initial value
spn_ds = tk.Spinbox(
    root, textvariable=var_ds, format="%.2f", from_=0.1, to=4.0, increment=0.1,
    command=lambda: change_slit_distance(var_ds.get()), width=5
    )
spn_ds.pack(side='left')

# Draw animation
set_axis()
anim = animation.FuncAnimation(fig, update, interval=50)
root.mainloop()
