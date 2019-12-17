"""
Creating visualization of the rock model
By: Yohanes Nuwara
"""

import matplotlib.pyplot as plt
import numpy.random as rnd
from matplotlib.patches import Ellipse

def rockmod(compos1, compos2):
    totalpopul = 300
    popul1 = int(compos1 * totalpopul)
    popul2 = int(compos2 * totalpopul)

    alpha1 = 0.1; alpha2 = 0.5
    width = 0.2
    height1 = width/alpha1; height2 = width/alpha2

    ells1 = [Ellipse(xy=rnd.rand(2)*10, width=width, height=height1, angle=rnd.rand()*360)
        for i in range(popul1)]

    ells2 = [Ellipse(xy=rnd.rand(2)*10, width=width, height=height2, angle=rnd.rand()*360)
        for i in range(popul2)]

    fig = plt.figure(0)
    fig = fig.add_subplot(111, aspect='equal')

    for e in ells1:
        def_ells1 = fig.add_artist(e)
        def_ells1 = e.set_clip_box(fig.bbox)
        def_ells1 = e.set_alpha(rnd.rand())
        def_ells1 = e.set_facecolor(rnd.rand(3))

    for e in ells2:
        def_ells2 = fig.add_artist(e)
        def_ells2 = e.set_clip_box(fig.bbox)
        def_ells2 = e.set_alpha(rnd.rand())
        def_ells2 = e.set_facecolor(rnd.rand(3))

        plot = fig.set_xlim(0, 10)
        plot = fig.set_ylim(0, 10)

    return(def_ells1, def_ells2, fig, plot)

    # plt.suptitle('Rock Elastic Property Model', fontsize=14, fontweight='bold')
    # plt.text(-0.1, -0.1, 'matplotlib', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
    #
    # plt.show()
    # return()
