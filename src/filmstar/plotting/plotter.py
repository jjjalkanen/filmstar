"""Interactive 3D plotter for mesh spline surfaces."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider


def update_plot(projection, tau):
    """Draw all faces; alpha rotates the mesh in the horizontal plane."""
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    fig.subplots_adjust(bottom=0.20)
    alpha_ax = fig.add_axes([0.15, 0.07, 0.7, 0.03])
    alpha_slider = Slider(alpha_ax, 'Alpha', 0, 2*np.pi, valinit=0)
    ax.set(xlabel='X', ylabel='Y', zlabel='Z', title=projection.title)
    ax.view_init(elev=28, azim=-55)
    artists = []
    colors = plt.get_cmap('tab10').colors

    def draw(alpha):
        for artist in artists:
            artist.remove()
        artists.clear()
        surfaces = projection.surface_func(tau, alpha)
        for k, surface in enumerate(surfaces):
            xx, yy, zz = np.moveaxis(surface, -1, 0)
            artists.append(ax.plot_surface(xx, yy, zz, color=colors[k % len(colors)],
                                           rcount=len(xx), ccount=xx.shape[1],
                                           linewidth=0, antialiased=True, alpha=0.96))
            for boundary in (surface[0], surface[-1], surface[:, 0], surface[:, -1]):
                artists.extend(ax.plot(*boundary.T, color='#333333', lw=0.65))
        return np.concatenate([part.reshape(-1, 3) for part in surfaces])

    points = draw(alpha_slider.val)
    radius = np.max(np.linalg.norm(points[:, :2], axis=1))*1.05
    zlo, zhi = np.min(points[:, 2]), np.max(points[:, 2])
    pad = max((zhi-zlo)*0.06, 0.1)
    ax.set(xlim=(-radius, radius), ylim=(-radius, radius), zlim=(zlo-pad, zhi+pad))
    ax.set_box_aspect((1, 1, 0.8))

    def update(alpha):
        draw(alpha)
        fig.canvas.draw_idle()

    alpha_slider.on_changed(update)
    # Keep the widget alive as long as the figure is open.
    fig._spline_slider = alpha_slider
    fig.canvas.draw()
    return fig
