import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


def fig_save_many(f, name, types=[".png",".pdf"], dpi=200):
    for ext in types:
        f.savefig(name+ext, dpi=dpi)
    return

def paper_single(TW = 6.64, AR = 0.74, FF = 1.):
    '''paper_single(TW = 6.64, AR = 0.74, FF = 1.)
    TW = 3.32
    AR = 0.74
    FF = 1.
    #mpl.rc('figure', figsize=(4.5,3.34), dpi=200)
    mpl.rc('figure', figsize=(FF*TW, FF*TW*AR), dpi=200)
    mpl.rc('figure.subplot', left=0.18, right=0.97, bottom=0.18, top=0.9)
    mpl.rc('lines', linewidth=1.0, markersize=4.0)
    mpl.rc('font', size=9.0, family="serif", serif="CM")
    mpl.rc('xtick', labelsize='small')
    mpl.rc('ytick', labelsize='small')
    mpl.rc('axes', linewidth=0.75)
    mpl.rc('legend', fontsize='small', numpoints=1, labelspacing=0.4, frameon=False) 
    mpl.rc('text', usetex=True) 
    mpl.rc('savefig', dpi=300)
    '''
    mpl.rc('figure', figsize=(FF*TW, FF*TW*AR), dpi=100)
    mpl.rc('figure.subplot', left=0.15, right=0.95, bottom=0.15, top=0.92)
    mpl.rc('lines', linewidth=1.75, markersize=8.0, markeredgewidth=0.75)
    mpl.rc('font', size=18.0, family="serif", serif="CM")
    mpl.rc('xtick', labelsize='small')
    mpl.rc('ytick', labelsize='small')
    mpl.rc('xtick.major', width=1.0, size=8)
    mpl.rc('ytick.major', width=1.0, size=8)
    mpl.rc('xtick.minor', width=1.0, size=4)
    mpl.rc('ytick.minor', width=1.0, size=4)
    mpl.rc('axes', linewidth=1.5)
    mpl.rc('legend', fontsize='small', numpoints=1, labelspacing=0.4, frameon=False) 
    mpl.rc('text', usetex=True) 
    mpl.rc('savefig', dpi=300)
    
    
    
def paper_single_mult_ax(nrows=1, ncols=1, **kwargs):
    #import matplotlib as mpl
    paper_single(FF=max(nrows,ncols))
    f, ax = plt.subplots(nrows=nrows, ncols=ncols, **kwargs)
    plt.minorticks_on()
    ylocator6 = plt.MaxNLocator(5)
    xlocator6 = plt.MaxNLocator(6)
    if len(ax.shape) > 1:
        for axrow in ax:
            for axcol in axrow:
                axcol.xaxis.set_major_locator(xlocator6)
                axcol.yaxis.set_major_locator(ylocator6)
    else:
        for axcol in ax:
            axcol.xaxis.set_major_locator(xlocator6)
            axcol.yaxis.set_major_locator(ylocator6)
    return f, ax
    
    
def paper_single_ax(TW = 6.64, AR = 0.74, FF = 1.):
    #import matplotlib as mpl
    paper_single(TW=TW, AR=AR, FF=FF)
    f = plt.figure()
    ax = plt.subplot(111)
    plt.minorticks_on()
    ylocator6 = plt.MaxNLocator(5)
    xlocator6 = plt.MaxNLocator(6)
    ax.xaxis.set_major_locator(xlocator6)
    ax.yaxis.set_major_locator(ylocator6)
    return f, ax

def paper_double_ax():
    #import matplotlib as mpl
    paper_single(TW = 12)
    f = plt.figure()
    ax = plt.subplot(111)
    plt.minorticks_on()
    ylocator6 = plt.MaxNLocator(5)
    xlocator6 = plt.MaxNLocator(6)
    ax.xaxis.set_major_locator(xlocator6)
    ax.yaxis.set_major_locator(ylocator6)
    return f, ax

def paper_double_mult_ax(nrows=1, ncols=1, setticks=True, **kwargs):
    #import matplotlib as mpl
    paper_single()
    
    
    TW = 6.97*2
    AR = 0.74
    FF = 1.
    mpl.rc('figure', figsize=(FF*TW, FF*TW*AR), dpi=200)
    mpl.rc('figure.subplot', left=0.1, right=0.97, bottom=0.1, top=0.97)
    mpl.rc('font', size=24.0, family="serif", serif="CM")
    
    f, ax = plt.subplots(nrows=nrows, ncols=ncols, **kwargs)
    plt.minorticks_on()
    if setticks:
        ylocator6 = plt.MaxNLocator(5)
        xlocator6 = plt.MaxNLocator(6)
        if len(ax.shape) > 1:
            for axrow in ax:
                for axcol in axrow:
                    axcol.xaxis.set_major_locator(xlocator6)
                    axcol.yaxis.set_major_locator(ylocator6)
        else:
            for axcol in ax:
                axcol.xaxis.set_major_locator(xlocator6)
                axcol.yaxis.set_major_locator(ylocator6)
    return f, ax




def hist2d(ax, xdat, ydat, xyrange, bins=[100,100], thresh=2, cmap=plt.cm.Greys, log=False, scatterother=False):
    import scipy

    tt = ax.get_aspect()

    # histogram the data
    hh, locx, locy = scipy.histogram2d(xdat, ydat, range=xyrange, bins=bins)
    mhh = np.mean(hh)
    shh = np.std(hh)
    if log:
        lhh = np.log10(hh)
    else:
        lhh = hh
    posx = np.digitize(xdat, locx)
    posy = np.digitize(ydat, locy)


    #select points within the histogram
    ind = (posx > 0) & (posx <= bins[0]) & (posy > 0) & (posy <= bins[1])
    hhsub = hh[posx[ind] - 1, posy[ind] - 1] # values of the histogram where the points are
    xdat1 = xdat[ind][hhsub < thresh] # low density points
    ydat1 = ydat[ind][hhsub < thresh]
    lhh[hh  < thresh] = np.nan # fill the areas with low density by NaNs

    ar = (0.6/0.65)*(np.diff(xyrange[0])/np.diff(xyrange[1]))[0]
    c = ax.imshow(np.flipud(lhh.T),extent=np.array(xyrange).flatten(), interpolation='none', cmap=cmap, aspect=ar)  
    
    ax.set_aspect(tt)
    
    if scatterother:
        ax.plot(xdat1, ydat1, 'k,')    
    
    
    return c


def make_ax3():
    paper_single(TW=8, AR=0.9)
    f = plt.figure()
    
    from matplotlib.ticker import NullFormatter, MaxNLocator

    nullfmt   = NullFormatter()         # no labels

    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.6
    bottom_h = bottom+height+0.02
    left_h = left+width+0.02

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    ax = plt.axes(rect_scatter)
    plt.minorticks_on()
    axx = plt.axes(rect_histx)
    plt.minorticks_on()
    axy = plt.axes(rect_histy)
    plt.minorticks_on()

    # no labels
    axx.xaxis.set_major_formatter(nullfmt)
    axy.yaxis.set_major_formatter(nullfmt)
    
    
    axy.xaxis.set_major_locator(MaxNLocator(3))
    axx.yaxis.set_major_locator(MaxNLocator(3))
    
    return f,ax,axx,axy



def nanhist(x,**kwargs):
    
    n,b,p = plt.hist(x[np.isfinite(x)],**kwargs)
    return n,b,p



