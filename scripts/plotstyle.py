import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

def set_style():

    # set font size ~appropriate to
    # a figure going in a paper
    sns.set_context("paper")
    
    # set font to be serif, not sans
    sns.set(font='serif')
    
    # white background, specify the
    # font family
    sns.set_style("ticks", 
       {"font.family": "serif"})#,
        #"font.serif": ["Times", "Palatino", "serif"]})

    # no axis frame
    plt.axes(frameon=False)

    # 26, 20, 14, 14, 14, "serif", "stix"
    params = {'axes.labelsize':32, 'axes.titlesize': 26, 'font.size':24,
              'xtick.labelsize': 24, 'ytick.labelsize': 24, 'font.family':"serif", 'mathtext.fontset':'stix'}
    #mpl.rcParams.update(params)

    return params

def set_aspect(ax, aspect="square"):

    if aspect == 'square':
        ax.set_aspect(np.abs((np.diff(ax.get_xlim()))/(np.diff(ax.get_ylim()))))

    return
