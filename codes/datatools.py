# -*- coding: utf-8 -*-
import numpy as np, scipy, networkx as nx, scipy.linalg as la
from math import *
import bisect as bisct,random
import matplotlib.pyplot as plt
import cProfile,pstats,sys,time,itertools,re
import inspect
from copy import deepcopy

###======================== DATA ANALYSIS TOOLSET ======================###


profiler=cProfile.Profile()

def setfigpath(path):
    '''Set path for saving figures.'''
    import os
    import matplotlib as mpl
    mpl.rcParams["savefig.directory"] = os.path.dirname(path)


def code_debugger(skip=0):
    import code
    import inspect
    stack=inspect.stack()
    def trace():
        print('\n'.join([' '.join([str(x[y]) for y in [1, 2]]) + '\n      ' + str(not x[4] or x[4][0]).strip() for x in
                         stack if x[4]]))
    dic = {}
    dic.update(stack[1+skip][0].f_globals)
    dic.update(stack[1+skip][0].f_locals)
    dic['trace']=trace
    code.interact(local=dic)


def rint(nb):
    return int(round(nb))
    
def kth_diag_indices(a, k):
    rows, cols = np.diag_indices_from(a)
    if k < 0:
        return rows[-k:], cols[:k]
    elif k > 0:
        return rows[:-k], cols[k:]
    else:
        return rows, cols


class MutableInt(object):
    def __init__(self,val):
        self.val = val
    def __add__(self,i):
        return self.val+i
    def __eq__(self,i):
        return self.val==i


def auto_subplot(plt,nbpanels,panel,projection=None,rows=None,return_all=0):
    i=panel.val
    if rows is None:
        panels_per_row=np.ceil(np.sqrt(nbpanels) )
    else:
        panels_per_row=np.ceil(nbpanels/rows).astype('int')
    nbrows=np.ceil(nbpanels*1./panels_per_row)
    ax=plt.subplot(100*nbrows+panels_per_row*10+i +1 ,projection=projection)
    panel.val+=1
    if return_all:
        return ax,nbrows,panels_per_row
    return ax

def gaussian(x,mean,var):
    return np.exp(-(x-mean)**2/2/var)/np.sqrt(2*np.pi*var)

def lognormal(x,mean,var):
    return np.exp(-(np.log(x)-mean)**2/2/var)/np.sqrt(2*np.pi*var)/x

def ifelse(cond,res,els):
    if cond:
        return res
    return els

def prolog(fname,col=2):
    stats=[[i.code ,i.totaltime,i.inlinetime,i.callcount,i.reccallcount] for i in profiler.getstats()]
    stats=sorted(stats,key=lambda e:e[col],reverse=1)
    with open(fname,'w') as prolog:
        for i in stats:
            if not i:
                continue
            st=' '.join([str(z) for z in  i])
            prolog.write(st+'\n')


def ADD_FOLDER_TO_PATH(folder,relative=1,depth=1):
    import os, sys, inspect
    if relative:
        cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile(
            inspect.currentframe(depth) ))[0])) +folder
    else:
        cmd_folder = os.path.realpath(os.path.abspath(folder))
    if cmd_folder not in sys.path:
        sys.path.insert(0, cmd_folder)



class Path(str):
    '''Strings that represent filesystem paths.
    Overloads __add__:
     - when paths are added, gives a path
     - when a string is added, gives a string'''
    def __add__(self,x):
        import os
        if isinstance(x,Path):
            return Path(os.path.normpath(os.path.join(str(self),x)))
        return os.path.normpath(os.path.join(str(self),x))

    def norm(self):
        import os
        return Path(os.path.normpath(str(self)))

    def osnorm(self):
        """Deal with different separators between OSes."""
        import os
        if os.sep=='/' and "\\" in str(self):
            return Path(os.path.normpath(str(self).replace('\\','/' )))
        elif os.sep=='\\' and "/" in str(self):
            return Path(os.path.normpath(str(self).replace('/','\\' )))
        else:
            return self.norm()

    def prev(self):
        import os
        lst=self.split()
        path=os.path.join(lst[:-1])
        return path.osnorm()

    def split(self):
        """"""
        import os
        lst=[]
        cur=os.path.split(self.norm())
        while cur[-1]!='':
            lst.insert(0,cur[-1])
            cur=os.path.split(cur[0])
        return lst

    def mkdir(self,rmdir=False):
        """Make directories in path that don't exist. If rmdir, first clean up."""
        import os
        if rmdir:
            os.rmdir(str(self))
        cur=Path('./')
        for intdir in self.split():
            cur+=Path(intdir)
            if not os.path.isdir(cur):
                os.mkdir(cur)

    def copy(self):
        return Path(self)

    def strip(self):
        '''Return string without final / or \\ to suffix/modify it.'''
        return str(self).strip('\/')

def wplot(typ,*args,**kwargs):
    USE_MYHIST=kwargs.get('USE_MYHIST',1)
    if kwargs.pop('newfig',0):
        plt.figure()
    save=kwargs.pop('save',0)
    hold=kwargs.pop('hold',False)
    title=kwargs.pop('title',False)
    legs=kwargs.pop('legends',False)
    invert=kwargs.pop('invert',0)
    projection=kwargs.pop('projection','2d')
    if not legs:
        legs=kwargs.pop('legend',False)
    if not legs:
        legs=kwargs.pop('leg',False)
    lkw=kwargs.pop('legopt',{})
    if 'xlabel' in kwargs:
        plt.xlabel(kwargs.pop('xlabel'))
    if 'ylabel' in kwargs:
        plt.ylabel(kwargs.pop('ylabel'))
    if 'labels' in kwargs:
        labs=kwargs.pop('xlabel')
        plt.xlabel(labs[0])
        plt.ylabel(labs[1])

    ax=plt.gca()

    if 'log' in kwargs:
        lg=kwargs.pop('log')
        if not isinstance(lg,str):
            if typ!= 'hist':
                kwargs['log']='xy'
            else:
                kwargs['log']=lg
        else:
            if 'x' in lg:
                ax.set_xscale('log')
            if typ!= 'hist' or USE_MYHIST:
                if 'y' in lg :
                    ax.set_yscale('log')
                if 'z' in lg:
                    ax.set_zscale('log')
            if typ=='hist':
                if 'y' in lg:
                    kwargs['log']=1
                if 'x' in lg:
                    if min(args[0])<=0:
                        args=([ar for ar in args[0] if ar>0],)+args[1:]
                    if not 'bins' in kwargs or isinstance(kwargs['bins'],int):
                        kwargs['bins']=np.logspace(log10(min(args[0])),
                            log10(max(args[0])),kwargs.get('bins',100),10)
    elif typ=='hist' and min(args[0])>0 and max(args[0])/min(args[0])>10**3:
        kwargs['log']=1
    if 'xs' in kwargs :
        xs=kwargs.pop('xs')
        args=list(itertools.chain.from_iterable([[xs[:len(a)],a] for a in args]))
    if not args:
        return
    handle =None
    force_zlim=kwargs.pop('force_zlim',False)

    while handle is None:
        try:
            if typ=='plot':
                handle=plt.plot(*args,**kwargs)
            if typ =='hist':
                if not 'bins' in kwargs:
                    kwargs['bins']=100
                if kwargs.pop('cleverbins',False):
                    if isinstance(kwargs['bins'],int):
                        kwargs['bins']=min(kwargs['bins'],len(args[0])/20)
                kwargs['hold']=1
                kwargs['density']=kwargs.pop('density',kwargs.pop('normed',1))
                if USE_MYHIST:
                    def myhist(a,**kws):
                        kw={}
                        plt.plot([],[])
                        kw.update(kws)
                        h=None
                        while h is None:
                            try:
                                h,b=np.histogram(a,**kw)
                            except TypeError as e:
                                # print e
                                del kw[str(e).strip().split(' ')[-1].replace("'",'').replace('"','') ]
                    #(b[1:]+b[:-1])/2,
                        handle=None

                        kw={}
                        kw.update(kws)
                        histtype=kw.get('histtype',None)

                        if kwargs.get('cumulative',0):
                            h=np.cumsum(h)
                            h/=h[-1]
                        while handle is None:
                            try:
                                if histtype=='step' and 0:
                                    x=b[:-1]
                                    y=h[:]
                                    if kw.get('log',0):
                                        x,y=x[x>0],y[x>0]
                                        x,y=x[y>0],y[y>0]
                                    handle=plt.plot( x, y, **kw )
                                else:
                                    handle=plt.bar(  b[:-1], h,b[1:]-b[:-1], **kw )
                            except AttributeError as e:
                                # print e
                                del kw[str(e).strip().split(' ')[-1].replace("'",'').replace('"','') ]
                        # plt.ylim(ymin=np.min(h[h>0]) )
                        # plt.autoscale(enable=1,axis='y')
                        return handle
                    handle=[myhist(a,**kwargs) for a in args ]
                else:
                    handle=[plt.hist(a,**kwargs)[2][0] for a in args]

            if typ=='scatter':
                nbaxes=2
                if projection=='3d':
                    nbaxes=3
                if len(args)>nbaxes:
                    handle=[]
                    for idx in range(len(args)/nbaxes):
                        handle.append(ax.scatter(*args[idx:idx+nbaxes],**kwargs))
                else:
                    handle=ax.scatter(*args,**kwargs)
                if force_zlim:
                    Z=args[2]
                    ax.set_zlim(bottom=min(Z), top=max(Z))

            if typ=='error':
                handle=plt.errorbar(*args,**kwargs)
            if typ=='errorfill':
                from mpltools import special as pltspecial
                handle=pltspecial.errorfill(*args,**kwargs)
            if typ=='contour':
                handle=plt.contour(*args,**kwargs)
            if typ=='contourf':
                handle=plt.contourf(*args,**kwargs)
            if typ=='bar':
                handle=plt.bar(*args,**kwargs)
        except AttributeError as e:
            #print e
            del kwargs[str(e).strip().split(' ')[-1] ]

        #except e:
            #print "Cannot plot:", sys.exc_info()[0],e
            #return

    if kwargs.get('colorbar',0):
        plt.colorbar(handle)
    if title:
        plt.title(title)
    if legs:
        plt.legend(handle,legs,**lkw)
    if invert:
        if 'x' in invert:
            plt.xlim(plt.xlim()[::-1])
        elif 'y' in invert:
            plt.ylim(plt.ylim()[::-1])
    if save:
        plt.savefig(save)
        plt.clf()
        plt.close()
        return
    elif not hold:
        plt.show()
    return handle

def errorbar(*args,**kwargs):
    return wplot('error',*args,**kwargs)
def errorfill(*args,**kwargs):
    return wplot('errorfill',*args,**kwargs)
def plot(*args,**kwargs):
    return wplot('plot',*args,**kwargs)
def hist(*args,**kwargs):
    return wplot('hist',*args,**kwargs)
def scatter(*args,**kwargs):
    return wplot('scatter',*args,**kwargs)
def scatter3d(*args,**kwargs):
    ax=plt.gcf().add_subplot(111,projection='3d')
    return wplot('scatter',*args,projection='3d',**kwargs)
def contour(*args,**kwargs):
    return wplot('contour',*args,**kwargs)
def contourf(*args,**kwargs):
    return wplot('contourf',*args,**kwargs)
def bar(*args,**kwargs):
    return wplot('bar',*args,**kwargs)

def cumulplot(*args,**kwargs):
    x=sorted(args[0])
    cx = np.linspace(0, 1, len(x))
    return wplot(kwargs.pop('plottype','plot'),x,cx,*args[1:],**kwargs)

def nonan(mat):
    return mat[np.isnan(mat)==False]
def nozero(mat):
    return mat[mat!=0]
def noinf(mat):
    return mat[np.logical_and(mat!=np.inf,mat!=-np.inf)]

