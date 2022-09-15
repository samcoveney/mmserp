# plotting function

import numpy as np
import pyvista as pv
pv.set_plot_theme("document")

from matplotlib.colors import LinearSegmentedColormap

def Tri_convert(Tri):
    """Creates a Tri representation suitable for pyvista."""
    return np.hstack([ np.full(Tri.shape[0], 3)[:,None] , Tri ])


def plot_fields(X, Tri, fields, titles, cmap = "jet", alpha = 1.0):
    """Plot field."""
    
    P = len(fields)

    # plot using pyvista
    plotter = pv.Plotter(shape = (1,P))

    # mesh and points
    plt_surf = pv.PolyData(X, Tri_convert(Tri))

    for i in range(P):

        plotter.subplot(0,i)
    
        clim = [fields[i].min(), fields[i].max()]
        scalars = fields[i]
        title = titles[i]

        # Controlling the text properties
        sargs = {"title": title}

        # plot mesh
        plotter.add_mesh(plt_surf.copy(), scalars = scalars, clim = clim, cmap = cmap, scalar_bar_args = sargs)

    plotter.link_views()
        
    return plotter


#        if i == 2:
#            clim = [CV_max.min(), CV_max.max()]
#            scalars = CV_max
#            title = "CV_max"


def plot_field_and_points(X, Tri, scalars = None, cmap = "jet", clim = None, title = "?", points = None, points_scalars = None, points_size = 25, fmt = "%.1f", alpha = 1.0, above_color = None):
    """Plot field on mesh and fie for paper, setup for Jupyter notebooks and png plots."""
    
    pv.rcParams['multi_rendering_splitting_position'] = 0.40
    
    # Controlling the text properties
    sargs = dict(
        title_font_size=int(20*2.75),
        label_font_size=int(16*2.75),
        shadow=True,
        n_labels=5,
        italic=False,
        fmt=fmt,
        font_family="arial",
        vertical = True,
        #position_x=0.75,
        #position_y=0.90,
        position_x=10, # will deliberatively hide scalar bar of first image
        position_y=0.05,
        height = 0.85,
        title = title + "_dummy")

    # plot using pyvista
    plotter = pv.Plotter(shape = "1|1", window_size = (1850,700), border = False)

    # mesh and points
    plt_surf = pv.PolyData(X, Tri_convert(Tri))
    
    if points is not None:
        plt_points = pv.PolyData(points)

    if (clim is None) and (scalars is not None):
        clim = [scalars.min(), scalars.max()]

    # loop over views of the mesh
    for i in range(2):

#        # create colorbar for second plot only, if scalars defined
#        if i == 1 and (scalars is not None):
#            stitle = title
#        else:
#            stitle = None
        
        plotter.subplot(i)

        if i == 1:
            sargs["title"] = title
            sargs["position_x"] = 0.75 

        # plot mesh
        if (scalars is not None):
            plotter.add_mesh(plt_surf.copy(), scalars = scalars,\
                             clim = clim, cmap = cmap, scalar_bar_args = sargs, above_color = above_color)
        else:
            plotter.add_mesh(plt_surf, color = "white", opacity = alpha)

        
        # plot points
        if points is not None:
            if (points_scalars is not None):
                plotter.add_points(plt_points.copy(), render_points_as_spheres = True, scalars = points_scalars,\
                               clim = clim, cmap = cmap, point_size = points_size, scalar_bar_args = sargs)
            else:
                plotter.add_points(plt_points, render_points_as_spheres = True, point_size = points_size, color = "black")

        # setup view -----------------------------
        if i == 0: plotter.view_xz()
        if i == 1: plotter.view_xz(negative = True)

        plotter.camera.Zoom(1.58)

        # offset second plot after rotation
        if i > 0:
            tmp = plotter.camera_position
            val = list(tmp[1])
            val[0] = val[0] - 17500
            tmp_new = [tmp[0], tuple(val), tmp[2]]
            plotter.camera_position = tmp_new

    return plotter


def colormap(cmap, b, t):
    """Can be used to truncate an existing colormap"""
    n = 500
    cb   = np.linspace(b, t, n)
    cm = cmap( cb )
    #cm[0:50,3] = np.linspace(0.80,1.0,50) # adjust alphas (transparencies)
    new_cmap = LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=b, b=t), cm )
    #new_cmap.set_under(color=myGrey())    
    new_cmap.set_over(color="black")    
    return new_cmap



