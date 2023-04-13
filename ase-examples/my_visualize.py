#!/usr/bin/env python
# coding: utf-8

# In[ ]:


def print_xyz(slab):
    import matplotlib.pyplot as plt
    from ase.visualize.plot import plot_atoms

    print(slab)
    fig, axarr = plt.subplots(1, 4, figsize=(15, 5))
    
    plot_atoms(slab, axarr[0], radii=0.3, rotation=('0x,0y,0z'))
    plot_atoms(slab, axarr[1], radii=0.3, rotation=('-90x,0y,0z'))
    plot_atoms(slab, axarr[2], radii=0.3, rotation=('-90x,-90y,0z'))
    plot_atoms(slab, axarr[3], radii=0.3, rotation=('-45x,0y,0z'))
    
    axarr[0].set_title("Top view")
    axarr[0].set_xlabel("X-axis, [$\mathrm{\AA}$]")
    axarr[0].set_ylabel("Y-axis, [$\mathrm{\AA}$]")
    
    axarr[1].set_title("Side view")
    axarr[1].set_xlabel("X-axis, [$\mathrm{\AA}$]")
    axarr[1].set_ylabel("Z-axis, [$\mathrm{\AA}$]")
    
    axarr[2].set_title("Side view")
    axarr[2].set_xlabel("Y-axis, [$\mathrm{\AA}$]")
    axarr[2].set_ylabel("Z-axis, [$\mathrm{\AA}$]")
    
    axarr[3].set_title("Diagonal view")
    axarr[3].set_axis_off()
    
    return

