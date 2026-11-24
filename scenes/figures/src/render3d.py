import numpy as np, pyvista as pv
from matplotlib.colors import LinearSegmentedColormap
pv.OFF_SCREEN = True
from style import SCENES, FIGURES
R = str(SCENES / 'review-export') + '/'
OUT = str(FIGURES) + '/'
ramp = ['#f4f8fd', '#cde2fb', '#9ec5f4', '#6da7ec', '#3987e5', '#256abf', '#184f95', '#0d366b']
cmap = LinearSegmentedColormap.from_list('blue', ramp)
INK = '#0b0b0b'
SB = dict(title='wind speed (m/s)\n', title_font_size=30, label_font_size=28, n_labels=6, fmt='%.0f',
          color=INK, vertical=True, position_x=0.9, position_y=0.12, height=0.6, width=0.04)

def label(p, pt, text):
    p.add_point_labels([pt], [text], font_size=32, text_color=INK, shape_color='white', shape_opacity=0.85,
                       show_points=False, always_visible=True)

# ---- coupled scene
m = pv.read(R + 'joint-baseline.vtu')
m.cell_data['speed'] = np.linalg.norm(m.cell_data['U'], axis=1)
mp = m.cell_data_to_point_data()
hs = mp.clip_box((-230, 230, -470, 230, 0, 130), invert=False).slice(normal='z', origin=(0, 0, 32.1))
vs = mp.clip_box((-1, 1, -470, -150, 0, 90), invert=False).slice(normal='x', origin=(0, 0, 0))
b = pv.read(R + 'buildings.stl')
ground = pv.Plane(center=(0, -120, -0.05), direction=(0, 0, 1), i_size=460, j_size=700)
p = pv.Plotter(off_screen=True, window_size=(1920, 1080)); p.set_background('white')
p.add_mesh(ground, color='#ecebe7', lighting=False)
p.add_mesh(b, color='#a9a7a0')
p.add_mesh(hs, scalars='speed', cmap=cmap, clim=(0, 10), opacity=0.8, lighting=False, scalar_bar_args=SB)
p.add_mesh(vs, scalars='speed', cmap=cmap, clim=(0, 10), lighting=False, show_scalar_bar=False)
p.add_mesh(pv.Disc(center=(0, -280, 32.1), inner=11.5, outer=14.5, normal=(0, 1, 0), c_res=120), color='#eb6834', lighting=False)
label(p, (45, -290, 8), 'rotor (actuator disk), hub 32.1 m')
label(p, (0, -470, 110), 'vertical section through the wake (x = 0)')
label(p, (-200, 200, 45), 'horizontal slice at 32.1 m')
p.add_arrows(np.array([[60, 300, 40]]), np.array([[0, -1, 0]]), mag=70, color=INK)
label(p, (60, 330, 75), 'wind from N')
p.camera_position = [(620, -780, 430), (0, -140, 15), (0, 0, 1)]
p.camera.zoom(1.0)
p.screenshot(OUT + 'fig01-coupled-scene-3d.png'); p.close(); print('coupled done')

