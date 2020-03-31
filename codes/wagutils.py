#from imports import *
from datatools import *
from scipy.stats import linregress
import matplotlib.figure as mpfig
from scipy.optimize import curve_fit
import numpy as np,sys, pandas as pd, os, ast
import matplotlib as mpl, matplotlib.pyplot as plt
import matplotlib.figure as mpfig
figsize=np.array(mpfig.rcParams['figure.figsize'])


def smoothness(mat,bs=None,mask=None):
    Slive=mat.shape[0]
    if mask is None:
        mask=np.ones(mat.shape)
    source=np.array([[(i,j) for j in range(Slive) if j!=i ] for i in range(Slive) ])
    def get(mat,x):
        return np.array([ mat[i,j] if  mask[i,j] else np.nan for i,j in x.reshape((-1,2))  ])
    #
    right=get(mat,source[:,1:])-get(mat,source[:,:-1])
    bot=get(mat.T,source[:,1:])-get(mat.T,source[:,:-1])
    if bs is None:
        bs=np.std(offdiag(mat))
    x=list(right)+list(bot)
    xnan=np.isnan(x)
    return   np.mean(( np.abs(np.array(x)[~xnan])<bs/Slive**.5) )



def getranks(x,**kwargs):
    from scipy.stats import rankdata
    if 'axis' in kwargs:
        # code_debugger()
        return np.apply_along_axis(getranks,kwargs['axis'],x)
    return rankdata(x,kwargs.get('method','average'))
    # ranks= np.argsort(np.argsort(x,**kwargs),**kwargs)
    # return ranks



def get_cmap(i,N):
    import matplotlib.cm
    return matplotlib.cm.get_cmap('jet')( (i*1./N ) )

def offdiag(x):
    return x[np.eye(x.shape[0])==0]

def genpool(S,mu=0.6,sigma=0.2,gamma=0.3,sigK=.5,muK=3.5,sigr=0.3,mur=-.5,**kwargs):
    Ks=np.random.lognormal(muK,sigK,S)
    rs=np.random.lognormal(mur,sigr,S)
    mat = np.random.normal(0, 1, (S, S))
    if gamma==0:
        asymm=0
    else:
        asymm=(1- np.sqrt(1-gamma**2) )/gamma
    mat= (mat+asymm * mat.T) / np.sqrt(1+ asymm**2)
    mat=mu+sigma*(mat-np.mean(offdiag(mat) ))/np.std( offdiag(mat))
    np.fill_diagonal(mat,0)
    beta=-mat
    return rs,Ks, beta

def dosimu(Aij,Ks,rs=None,noise=1,evol=1,tmax=1000.,ts=None,**kwargs):
    if ts is None:
        ts= np.linspace(1,tmax,20)
    Aij=Aij.copy()
    x0=np.abs(np.random.normal(1.4,.4)*Ks/len(Ks))
    from scipy.integrate import odeint

    if rs is None:
        assert Aij[0,0]<0
        rs=-np.diag(Aij)*Ks
    np.fill_diagonal(Aij,0)
    
    def LV(x,t):
        return rs *x*(1 - x/Ks )+ x*np.dot(Aij,x) 
    result = odeint(LV,x0,ts)
    result*=(1+noise/10.*np.random.gamma(1,1,result.shape))
    result[result<10**-6]=0
    return ts,result


def find_eqs(beta,uninvadable=1,stable=1,largest=0):
    """Find equilibria. Here by convention, beta is negative, corresponds to K=1, and includes the diagonal -1"""
    from scipy.linalg import inv as lainv,eigvals,LinAlgError
    S=beta.shape[0]
    beta=beta.copy()
    assert np.max(np.abs(beta[0,0]+1))<10**-5
    pocket=[range(S)]
    eqs=[]
    while pocket:
        cur=list(pocket.pop(0))
        if not len(cur):
            break
        s=len(cur)
        betacur=beta[np.ix_(cur,cur)]
        try:
            eq=-np.dot(lainv(betacur),np.ones(s))
        except LinAlgError as e:
            continue
        res=np.zeros(S)
        res[cur]=eq
        notstab = stable and np.max(eigvals(betacur)) > 10 ** -10
        if s < S:
            notuninv = uninvadable and np.max(1 + np.dot(beta, res)[res == 0]) > 0
        else:
            notuninv = False
        feasible =not (eq<0).any()
        if not feasible or  notstab or notuninv:
            added=[z for z in [tuple(x for x in cur if x != cur[i]) for i in range(s)] if not z in pocket]
            pocket += added
        else:
            eqs.append(res)
            if largest:
                return eqs
    return eqs

def nonan(mat):
    return mat[np.logical_or(np.isnan(mat),np.isinf(mat))==False]


def cleverfit(xs,ys,year0=0,minr=0.25,clever=1,typ='mono',curve=0,debug=0):
    """Utility function to fit trajectories with either exponential or logistic growth"""
    def f(year0,x,a,b,c):
        return a + b*np.exp(-c*1.*np.abs(x-year0 ) )
    def monof(year0,x,a,b,c):
        return a / (1 - b * np.exp(-c * (x - year0)))
    fnc=None
    if typ=='mono':
        fnc=monof
    if typ=='exp':
        fnc=f
    def monofit(xs,ys,minr=0.25,maxr=3,year0=0):
        bprm=None
        if typ=='mono':
            bprm=-.3,.3,.99  #min(0,1-xs[-1]/max(xs[0],0.01) )
        if typ=='exp':
            bprm=0,ys[xs==year0].mean(),np.max(ys)
        p0=(np.median(ys), bprm[1], np.atleast_1d( (minr*9+maxr)/10.)[0] )
        return curve_fit(lambda *args: fnc(year0,*args),xs,ys,p0=p0,
                        bounds=((0,bprm[0],max(0,min(maxr,minr))),(np.max(ys),bprm[-1],maxr)) )
    tmp=[]
    def do(minr):
        p,var=monofit(xs,ys,minr,year0=year0)
        # relvar=np.diag(var)/p**2
        relvar=np.abs( var/np.multiply.outer(p,p) )
        if debug:
            tmp.append( (np.sum(relvar),np.sum(np.abs( (var/np.multiply.outer(p,p) ) [np.eye(len(p))==0 ]))))
        return p,var,np.sum(relvar)
    if len(set(xs))>1 and len(set(ys))>3: #need more points than parameters
        p, var, cost=do(minr)
    else:
        p=np.zeros(3)
        p[0]=np.mean(ys)
        var=np.zeros((3,3))
        if curve:
            return p, var ,pd.Series(data=ys,index=xs),np.arange(len(set(xs)))
        return p, var
    if p[-1]-minr <10**-3 and clever:
        import scipy.optimize as sopt
        try:
            res=sopt.minimize(lambda m:do(m)[-1],x0=minr,bounds=np.array([(0,1)]) )
        except:
            return cleverfit(xs,ys,year0=year0,minr=minr,clever=0,typ=typ,curve=curve,debug=debug)
        minr=res.x

    p, var, cost=do(minr)
    if debug:
        print var,var/np.multiply.outer(p,p)
        plt.figure()
        xx, yy = zip(*tmp)
        scatter(xx, yy, log='xy')
    if clever and len(set(xs))>2:
        #TEST IF GOING EXPONENTIALLY TO EXTINCTION
        pys=fnc(year0,xs,*p )
        half=(xs>np.min(xs))
        pexp,varexp=curve_fit(lambda *args: f(year0, args[0],0,args[1],args[2]), xs[half], ys[half], p0=(np.max(ys)/2,.5),
                      bounds=((0,min(0.4,1./len(xs)) ), ( np.max(ys), 10)))
        expys = f(year0, xs[half], 0, *pexp)
        def ddx(x):
            return np.log(x[1:]/x[:-1])
        def dist(x,y):
            return np.mean( (x-y)**2  )
            # return 1-np.corrcoef(ddx(x),ddx(y))[0,1]
        # print dist(pys[half],ys[half]),dist(expys,ys[half])
        if  dist(pys[half],ys[half])>2*dist(expys,ys[half]):
            p,var= np.array((0.0*np.min(ys),pexp[0],pexp[1])),varexp
            typ='exp'
            fnc=f
            # test=fnc(year0,xs,*p)
            # scatter(xs,ys,hold=1,log='y')
            # plot(pys,test,xs=xs,hold=1)
            # plt.legend(['Normal','test'])
            # plt.show()
            # code_debugger()


    if curve:
        pys=fnc(year0,xs,*p )
        pys=pd.Series(data=pys,index=xs)
        pys=pys.groupby(pys.index).first()
        eqyears=np.where( np.abs(pys/p[0]-1) <0.2  )[0]
        return p, var, pys, eqyears
    return p,var

def linfit(xs,ys,log=''):
    """Convenience function for linear regression and subsequent plot."""
    x,y=xs,ys
    if 'x' in log:
        x=np.log10(x)
    if 'y' in log:
        y = np.log10(y)
    slope, intercept, r, p, stderr = linregress(x,y)
    pxs=np.sort(xs)
    pys=np.sort(x)
    pys=slope*pys+intercept
    if 'y' in log:
        pys=10**(pys)
    desc=r'${:.2f} + {:.2f} x$ ($R^2$={:.2f}, p={:.2f})'.format(intercept,slope,r**2,p)
    if 'x' in log and 'y' in log:
        desc=r"${{{: .2f}}} \; x^{{{: .2f}}}$ ($R^2$={:.2f}, p={:.2g})".format(10**intercept,slope,r**2,p)
    return pxs,pys,desc,slope,intercept,r,p,stderr

def ks_nd(x,y, use=['sum', 'max'][1],log=0,debug=0):
    """N-dimensional Kolmogorov Smirnov-like test of x (data) against y (proposed distribution)"""
    bx,by=np.atleast_2d(x),np.atleast_2d(y)
    if bx.shape!=x.shape:
        bx,by=bx.T,by.T
    if log:
        for b in bx,by:
            b[b<0]=0
            b[b>0] = np.log(b[b > 0])
    rx=getranks(bx,axis=0)
    rjoint=getranks(np.concatenate([bx,by] ),axis=0   )
    rxy=rjoint[:rx.shape[0] ]
    Fx=rx*1./rx.shape[0]
    Fxy=rxy*1./rjoint.shape[0]
    if debug:
        code_debugger()
    if use =='max':
        return np.max(np.sum(np.abs(Fx-Fxy),axis=1))
    return np.mean(np.sum((Fx-Fxy)**2,axis=1)**0.5 )


def make_cc(S=8):
    """Make a beta matrix for CC Tradeoff structure"""
    bb = np.eye(S)
    m = 1.
    c = 2.
    x =np.array([ 1.])
    cs=[]
    for kk in range(1,S):
        if len(cs):
            c=cs[-1]
        ci = m / (1 - np.sum(x * (c - m) / c)) ** 2 + 2
        cs.append(ci)
        bb[kk, :kk] = [  (ci + cs[j]) / cs[j]* (cs[j] - m) / (ci - m)  for j in range(kk)]
        x = np.concatenate([x,[1 - np.sum(bb[kk, :kk] * x)]])
    return bb


def make_cctradeoff(eta,recursing=0,xprev=None,ftol=0.01,trial=0):
    """Make a beta (rescaled interaction) matrix according to CC Tradeoff structure that gives the desired etas (N/K)"""
    import scipy.linalg as la
    trueeta=eta.copy()
    eta=np.sort(eta)[::-1]
    # eta[0]=1
    eta/=np.max(eta)
    S=len(eta)
    def beta(c,m):
        b=np.zeros((S,S))
        i,j=np.indices(b.shape)
        b[j<i]=( np.add.outer(c,c)* np.multiply.outer( 1./m, m/c )  )[j<i]  #m  here is actually c-m
        b[np.isnan(b)]=0
        b[np.isinf(b)]=0
        np.fill_diagonal(b,1)
        return b
    def cost(x):
        c=x[:S]
        # m=np.clip(x[S:],0.0001,1)*c
        m = np.clip(np.exp(-x[S:]), 0, 1) * c
        # cost=(np.dot(la.inv(beta(c,m) ),np.ones(S)) -eta)
        cost= np.dot(beta(c,m),eta)-1
        # print c,m,la.norm(cost)*(1+np.sum(x[x<0])**2)
        # print 'Cost',la.norm(cost), c,m
        return la.norm(cost)*(1+np.sum(x[x<0])**2)
    import scipy.optimize as sopt
    success=0
    def geteta(c,m):
        b=beta(c,m)
        et=np.dot(la.inv(b),np.ones(S))
        return et
    tries=0
    if not recursing:
        xprev = np.linspace(1, S, S) ** 1.3 + np.random.random(S)
        xprev = np.concatenate([xprev, np.ones(S) * .1])
        for size in range(min(len(eta),4), len(eta)+1):
            indices=range(np.ceil(size/2.).astype('int') )+range(S-np.floor(size/2).astype('int'),S )
            indices2=[z+S for z in indices]
            print 'Recursing', size,# indices,len(indices)
            res=make_cctradeoff(eta[indices],recursing=1,xprev=xprev[indices+indices2],trial=trial )
            if res is None:
                print 'cc_tradeoff: FAIL!'
                return None
            c, m, b=res
            xprev[indices]= c
            xprev[indices2]=-np.log(m/c)
    else:
        while not success:
            x0 = np.linspace(1, S, S) ** 1.3 + np.random.random(S)
            x0 = np.concatenate([x0, np.ones(S) * .1])
            if not xprev is None:
                x0[:len(xprev)/2] = xprev[:len(xprev)/2] * np.random.normal(1,.1,len(xprev)/2)
                x0[S:S+len(xprev)/2] = xprev[len(xprev)/2:] * np.random.normal(1,.1,len(xprev)/2)
            res=sopt.minimize(cost,x0 ,
                              bounds=np.array((np.ones(2*S)*0.0001,np.concatenate([np.ones(S)*10**3,np.ones(S)*2*S] )) ).T )
            x = res.x
            c = x[:S]
            m = np.clip(np.exp(-x[S:]), 0, 1) * c
            print cost(x)#, np.dot(beta(c,m),eta)
            success= cost(x)<ftol
            tries+=1
            if tries>2:
                return None
        b=beta(c,m)
    idxs=np.argsort(np.argsort(trueeta)[::-1])
    c,m,b=c[idxs],m[idxs],b[np.ix_(idxs,idxs)]
    return c,m,b


def snap_beta(beta,eta,val=1,diag=1):
    """Return beta snapped to beta.eta=val"""
    betatmp=beta.copy()
    np.fill_diagonal(betatmp,0)
    lam= (val-diag*eta-np.dot(betatmp,eta))/(np.sum(eta**2) - eta**2 )
    beta_sn= beta+np.multiply.outer(lam,eta)
    np.fill_diagonal(beta_sn,diag)
    try:
        assert np.allclose(np.dot(beta_sn,eta),val )
    except:
        print '    snap_beta PROBLEM',np.dot(beta_sn,eta)
    return beta_sn


def make_worst(eta,bm,bs,trials=500000,invtemp=10. ):
    """Make matrix that is as far removed as possible from betatheo"""

    bm, bv, betatheo=hebbian_getbstats(np.ones((eta.shape[0],eta.shape[0])),eta,bm=bm,bv=bs**2)

    def bilin(etai,etaj,mat):
        def bi(x,A,B,C,D):
            return A*x[0]*x[1] + B*x[0] + C*x[1] +D
        res=curve_fit(bi,np.array([etai,etaj]),mat )
        # print res
        return res[0][0]

    def collision(phi,v1x,v2x,v1y,v2y):
        a=np.tan(phi)
        d2=2*(v1x-v2x + a*(v1y -v2y ))/(1+a**2)/2.
        pv1x=v1x-d2
        pv2x=v2x+d2
        pv1y=v1y-a*d2
        pv2y=v2y+a*d2
        return pv1x,pv2x,pv1y,pv2y

    beta=create_synthetic_network_LV_gam0(eta, bm , bs)
    off=offdiag(beta)
    offtheo=offdiag(betatheo)
    ii=np.multiply.outer( np.arange(len(eta)),np.ones(len(eta)) ).astype('int')
    i,j=offdiag(ii),offdiag(ii.T)
    matcons=np.array([ (i==k)*(eta[list(j)]) for k in range(len(eta))])
    veccons=1-eta
    etai,etaj=eta[list(i)],eta[list(j)]

    slop=bilin(etai,etaj,offtheo)
    # code_debugger()
    def dist(x):
        if True or 'bilin' in sys.argv:
            s=bilin(etai,etaj,x)
            return s*np.sign(slop)
        return  - np.sum(  np.abs(x-offtheo) )

    if 0:
        oldbeta = beta.copy()
        offold = offdiag(oldbeta)
        rej1,rej2=0,0
        preverr=0
        for trial in range(trials):
            if trial%10000==0:
                print trial
            i, j,k,l = np.random.choice( range(len(off)),replace=False,size=4)
            phi=np.random.uniform(0,np.pi*2)
            prevbetas=off[[i,j,k,l]]
            newbetas=np.array(collision(phi,*prevbetas))

            off2=off.copy()
            off2[[i, j, k, l]] = newbetas

            if np.random.random() < np.exp(invtemp * np.sum( - (offtheo[[i,j,k,l]] - prevbetas )*( newbetas-prevbetas)  )    ):
                err=np.mean((np.dot(matcons,off2)-veccons)**2)**.5

                if np.random.random()<np.exp( invtemp*(-err+preverr)):
                    off=off2
                    preverr=err
                else:
                    rej2+=1
            else:
                rej1+=1
        beta[np.diag(np.ones(eta.shape[0]) )==0 ]=off
        code_debugger()

    from scipy.optimize import minimize,NonlinearConstraint, LinearConstraint
    empbm=np.mean(off)
    empbs=np.std(off)
    if 0:
        eq_bm = {'type': 'eq', 'fun': lambda x:  np.mean(x) - empbm ,
         # 'jac': lambda x: np.ones(x.shape)*1./len(x)
                 }
        eq_bs = {'type': 'eq', 'fun': lambda x:  np.var(x) - empbs**2 ,
         # 'jac': lambda x: (2*x - 2.*empbm)/len(x)
                 }
        eq_eta = {'type': 'eq', 'fun': lambda x:np.sum( (np.dot(matcons,x)-veccons)**2 ),}
        res = minimize(dist, off, method='SLSQP',constraints = [eq_bm,eq_bs,eq_eta], options = {'ftol': 1e-9, 'maxiter':5000, 'disp': True},)
    else:
        b=empbm*np.ones(len(off))
        eq_bm=LinearConstraint(np.ones((len(off),len(off)) ),lb=b,ub=b )
        b=veccons
        eq_eta=LinearConstraint(matcons,lb=b,ub=b)
        eq_bs=NonlinearConstraint(lambda x: np.var(x),0,empbs**2 )
        eq_bm=NonlinearConstraint(lambda x: np.mean(x)-empbm,0,0 )
        res = minimize(dist, off, method='trust-constr',constraints = [eq_bm,eq_bs,eq_eta], options = {'maxiter':5000, 'disp': False},)

    beta[np.diag(np.ones(eta.shape[0])) == 0] = res.x

    # compmats = [offdiag(create_synthetic_network_LV_gam0(eta, bm, bs)) for z in range(100)]
    # scatter([np.std(b) for b in compmats], [bilin(etai, etaj, b) for b in compmats])
    # code_debugger()
    if np.std(res.x)<empbs:
        corr=snap_beta(create_synthetic_network_LV_gam0(eta, 0, (empbs**2-np.var(res.x) )**.5,val=0),eta,val=0,diag=0)
        beta=beta+corr
    return beta


def hebbian_getpairs(beta,eta,remove_zeros=0):
    eta= np.array(eta)
    S = eta.shape[0]
    indices = [(i, j) for i in range(S) for j in range(S) if i != j]
    if len(eta.shape)==1:
        etapairs=np.array([(eta[i],eta[j]) for i,j in indices]).T
    else:
        tetapairs=np.array([(eta[i,0],eta[j,0]) for i,j in indices]).T
        indices = [(i, j, t) for t in range(eta.shape[1]) for i in range(S) for j in range(S) if i!=j  ]
        etat=eta.T.ravel()
        etapairs=np.array([(etat[t*S+i],etat[t*S+j]) for i,j,t in indices]).T
        assert (etapairs[:,:tetapairs.shape[1]]==tetapairs).all()
        # code_debugger()
    betapairs=np.array([beta[i[0],i[1]] for i in indices]).T
    if remove_zeros:
        wh=np.where(np.logical_and(betapairs!=0,np.logical_and(np.min(etapairs,axis=0)>0,np.abs(betapairs)<1.5)) )[0]
        return etapairs[:,wh],betapairs[wh]
    return etapairs,betapairs


def hebbian_gettriples(beta,eta,full=0,remove_zeros=0,repeats=()):
    """Full not necessary for correlations in"""
    if repeats and not hasattr(repeats[0],'__iter__'):
        repeats=repeats,
    eta= np.array(eta)
    S = eta.shape[0]
    def cond(i,j,k):
        return ( i!=j or (0,1) in repeats) and ( i!=k or (0,2) in repeats) and ( j!=k or (1,2) in repeats)
    indices = [(i, j, k) for i in range(S) for j in range(S) for k in range(S) if cond(i,j,k) ]
    if len(eta.shape)==1:
        etatriples=np.array([(eta[i],eta[j],eta[k]) for i,j,k in indices]).T
    else:
        # tetapairs=np.array([(eta[i,0],eta[j,0]) for i,j in indices]).T
        indices = [(i, j,k, t) for t in range(eta.shape[1]) for i in range(S) for j in range(S)  for k in range(S) if cond(i,j,k) ]
        etat=eta.T.ravel()
        etatriples=np.array([(etat[t*S+i],etat[t*S+j],etat[t*S+k]) for i,j,k,t in indices]).T
        # assert (etapairs[:,:tetapairs.shape[1]]==tetapairs).all()
        # code_debugger()
    if full:
        betatriples=np.array([(beta[i[0],i[1]],beta[i[0],i[2]],beta[i[1],i[0]],beta[i[2],i[0]] ) for i in indices]).T
    else:
        betatriples=np.array([(beta[i[0],i[1]],beta[i[0],i[2]] ) for i in indices]).T
    if remove_zeros:
        # BAD IDEA: REMOVE ZEROS -- THIS WAAS SUPPOSED TO SPEED UP PROCESSING, BUT IT CHANGES AVERAGES AND RESULTS!
        wh=np.where(np.logical_and(np.sum((betatriples!=0),axis=0)==betatriples.shape[0],np.min(etatriples,axis=0)>0) )[0]
        return etatriples[:,wh],betatriples[:,wh]
    return etatriples,betatriples

def hebbian_loglike(beta,bm,bs):
    return np.mean( -(offdiag(beta)-bm)**2 /(2*bs**2) +.5 )

def hebbian_stablebm(eta):
    S=len(eta)
    def avgsum(x):
        return np.mean(np.sum(x, axis=0))
    qsum2 = avgsum(eta**2)
    Q1=np.mean(1/(qsum2-eta**2) )
    Q2= np.mean(eta /(qsum2-eta**2)   )
    bm= (Q1-Q2)/( S*np.mean(eta)*Q1-Q2  )
    return bm

def hebbian_etastar(eta,bm):
    def avgsum(x):
        return np.mean(np.sum(x, axis=0))
    return (1 - bm*avgsum(eta))/(1-bm)

def hebbian_convert(eta,bm,bs,forward=1,debug=0):
    '''Go between pool stats and survivor stats'''
    def avgsum(x):
        return np.mean(np.sum(x, axis=0))
    bv=bs**2
    qsum2 = avgsum(eta**2)
    c1=np.mean(1/(qsum2-eta**2) )
    c2= np.mean(eta /(qsum2-eta**2)   )
    A1= (c2-c1)  *np.mean(eta)
    A2=1- np.mean(eta)*( c1*avgsum(eta) -c2   )
    C=(1- np.mean(eta**2)*c1)
    if np.isnan([A1,A2,C,c1,c2]).any():
        if debug:
            code_debugger()
        else:
            raise
    if forward:
        bmp=bm*A2  - A1
        bvp=bv*C
        return {'bm_pool': bm, 'bs_pool': bs, 'bm_surv': bmp, 'bs_surv': bvp**.5}
    else:
        '''Get pool stats from survivor stats'''
        bmp= (bm + A1)/A2
        bvp=bv/C
        return {'bm_pool': bmp, 'bs_pool': bvp**.5, 'bm_surv': bm, 'bs_surv': bs}

def hebbian_getbstats(beta,eta,**kwargs):
    """Get pool mean and variance from observed beta and eta"""
    S=eta.shape[0]
    eta=eta.astype('float')
    betaerr=kwargs.get('betaerr', np.ones(beta.shape))**2
    if  kwargs.get('bm',None) is None:
        if kwargs.get('use_MF_bm',0 ):
            # print "USING MEANFIELD"
            bm_surv=(1./np.mean(eta)-1)/S
        elif kwargs.get('get_bm_from_beta',0):
            # print 'bm FROM DATA'
            ob=offdiag(beta)
            bm_surv = np.average(ob[~np.isnan(ob)],weights=1./offdiag(betaerr)[~np.isnan(ob)])
        else:
            # print "NOT MEANFIELD"
            bm_surv=hebbian_stablebm(eta)
        try:
            bloc=hebbian_convert(eta,bm_surv,1,forward=0)
        except Exception as e:
            print e
        bm=bloc['bm_pool']
    else:
        bm=kwargs['bm']
    if len(eta.shape)>1:
        etam=np.mean(eta,axis=1)
    else:
        etam=eta
    etai,etaj,etaij=np.multiply.outer(etam,np.ones(S)),np.multiply.outer(np.ones(S),etam),np.multiply.outer(etam,etam)
    sum1,sum2=np.mean(np.sum(eta,axis=0)),np.mean(np.sum(eta**2,axis=0))
    betatheo=bm - etaij/(sum2-etai**2) + (1- bm*(sum1 - etai) )/(sum2-etai**2) *etaj
    np.fill_diagonal(betatheo,1)
    if not 'bv' in kwargs:
        # bs_surv=np.mean(  offdiag((beta-betatheo)**2) )**.5
        # bv_surv = np.mean(offdiag((beta-betatheo)**2/betaerr**2))/np.mean(offdiag(1./betaerr**2))
        ob = offdiag(beta)
        bv_surv=np.average((offdiag(beta-betatheo)**2)[~np.isnan(ob)], weights=1./offdiag(betaerr)[~np.isnan(ob)] )
        """IMPORTANT :"""
        # code_debugger()
        bs_surv=bv_surv**.5
        bloc=hebbian_convert(eta, bm,bs_surv,forward=0)
        bv=bloc['bs_pool']**2
    else:
        bv=kwargs['bv']
    return bm,bv,betatheo


def hebbian(beta,eta,with_gamma=0):
    eta= np.array(eta)
    S = eta.shape[0]
    etapairs,betapairs=hebbian_getpairs(beta,eta)
    if with_gamma:
        gamma=np.corrcoef(offdiag(beta),offdiag(beta.T))[0,1]
    else:
        gamma=0
    # print betapairs[:10], beta[zip(*indices) ]
    bm, bv, betatheo=hebbian_getbstats(beta,eta)
    bs=bv**.5
    # if bm<0 :
    #     print 'WARNING: PROBABLY WRONG SIGN FOR BETA'
    def hebb(x,A,B1,B2):
        return (bm-A*x[0]*x[1] + B1*x[0] + B2*x[1] )
    popt,pcov=scipy.optimize.curve_fit(hebb,etapairs,betapairs,p0=(0,0,0) )
    # code_debugger()
    # if len(eta.shape)>1:
    #     eta=np.mean(eta,axis=tuple(range(1,len(eta.shape))) )
    def avgsum(x):
        return np.mean(np.sum(x, axis=0))
    qsum=avgsum(eta)*(S-1)/S
    qsum2=avgsum(eta**2)*(S-1.)/S
    Btheo=(1-bm*qsum )/qsum2
    Atheo=(1+gamma)*(1.- bm + gamma*Btheo * qsum ) /qsum2  #Seemingly a problem!
    return {'A':popt[0],'Bgamma':popt[1],'B':popt[2],'covbetaAB':pcov, 'A/theo': popt[0]/Atheo,
            'B/theo':popt[2]/Btheo, 'Atheo':Atheo,'Btheo':Btheo, 'Bgammatheo':Btheo*np.corrcoef(offdiag(beta),offdiag(beta.T))[0,1],
            'loglike':hebbian_loglike(beta,bm,bs) }

def hebbian_full(beta,eta,with_gamma=0):
    eta = np.array(eta)
    S = eta.shape[0]
    bm, bv, betatheo=hebbian_getbstats(beta,eta)
    def avgsum(x):
        return np.mean(np.sum(x, axis=0))
    dic={'Ctheo':np.mean(1./(avgsum(eta**2) - eta**2 ) ),  }
    dic.update( hebbian(beta, eta, with_gamma=with_gamma))
    etatriples,betatriples=hebbian_gettriples(beta-betatheo,eta,repeats=(1,2))
    # code_debugger()

    #
    def Hcov(x,C):
        djk=np.zeros(x.shape[1])
        djk[x[1]==x[2]]=1
        return bv*(-C*x[1]*x[2] +  djk )

    def qcov(x, C):
        djk = np.zeros(x.shape[1])
        djk[x[1] == x[2]] = 1
        return bv * (C + djk)

    xy=etatriples
    z=np.product(betatriples,axis=0)
    popt,pcov=scipy.optimize.curve_fit(Hcov,xy,z,p0=(0,) )
    qopt,qcov=scipy.optimize.curve_fit(qcov,xy,z,p0=(0,) )
    # test=linregress(np.product(etatriples[1:],axis=0), -z )
    #
    dic['C']=popt[0]
    # dic['C2']=test[0]/bv
    dic['Cstd']=pcov[0,0]**.5
    dic['Crelstd']=pcov[0,0]**.5/np.abs(popt[0])
    dic['naiveC']=-np.mean(z)/np.mean(np.product(xy,axis=0) )/bv
    dic['correl']=qopt[0]
    dic['correltheo']=-dic['Ctheo']*np.mean(np.product(xy,axis=0) )
    # dic['naiveC']=-np.mean(z/np.product(xy,axis=0) )/bv
    return dic

def hebbian_slope(beta,eta,with_gamma=0,bins=8,plots='A',etavar=0,**kwargs):
    """WARNING: A obtained here (by slope of slopes) has much wider distribution that that obtained by bilinear fit"""
    etapairs,betapairs=hebbian_getpairs(beta,eta)
    if not hasattr(bins,'__iter__'):
        if kwargs.get('log',0):
            bins=np.logspace(np.log10(np.min(eta)+10**-3), np.log10(np.percentile(eta,90)), bins)
        elif kwargs.get('adapt_bins',0):
            from scipy.cluster.hierarchy import fclusterdata
            nbins=bins
            handled=False
            while not handled and nbins>2:
                cl=fclusterdata(eta,nbins,'maxclust')
                means=np.sort([np.mean(eta[cl==i]) for i in sorted(set(cl))])
                bins=[0]+list( (means[1:]+means[:-1])  /2 )+[np.max(eta)*1.01]
                handled=min(np.histogram(eta.ravel(),bins=bins)[0])>1
                nbins-=1
        else:
            bins=np.linspace(np.min(eta[eta>0])*.99,np.max(eta)*1.01, bins+1)
    # MEAN AND CORRELATION
    bm, bv, betatheo = hebbian_getbstats(beta, eta,**kwargs)
    oldbm=bm
    etatriples,betatriples=hebbian_gettriples(beta- betatheo,eta, repeats=(1,2),full=1)

    betaerr=kwargs.get('betaerr',np.ones(beta.shape))**2
    tmp,betaerrpairs=hebbian_getpairs(betaerr,eta)
    tmp,betaerrtriples=hebbian_gettriples(betaerr,eta, repeats=(1,2),full=1)

    ##
    idx,legends=0,[]
    table=[]
    # def fun(x,  B):
    #     return bm + B *  x
    for bin1,bin2 in zip(bins[:-1],bins[1:]):
        wh=np.where(np.logical_and(etapairs[0]>=bin1,etapairs[0]<=bin2 ))[0]
        twh=np.where(np.logical_and(etatriples[0]>=bin1,etatriples[0]<=bin2 ))[0]

        if  len(wh)<2:
            continue
        popt, pcov = scipy.optimize.curve_fit(lambda x,a,b:a*x+b, etapairs[1,wh], betapairs[wh], p0=(0,0 ),sigma=betaerrpairs[wh])
        if plots=='all':
            plt.scatter(etapairs[1,wh],betapairs[wh],alpha=.5,marker='x',color=get_cmap(idx,len(bins)-1))
            def rge(x):
                return np.min(x),np.max(x)
            plt.plot(rge(etapairs[1,wh]),bm+popt[0]*np.array(rge(etapairs[1,wh])),color=get_cmap(idx,len(bins)-1) )
            # rge = np.array([0, 1])
            # plt.plot(rge, popt[0] + popt[1] * rge,color=get_cmap(idx,len(bins)-1))
        idx+=1
        dic={'legend':'{:.2f}'.format((bin1+bin2)/2), 'slope':popt[0],'const':popt[1],'nonempty':np.sqrt(bin1*bin2),
             'binstart':bin1,'binend':bin2, 'stderrslopes':pcov[0,0]**.5 / np.sqrt(len(set(betapairs[wh]) )),
             'errslope':pcov[0,0]**.5, 'meanetai':np.mean(etapairs[0,wh]),'meanetai2':np.mean(etapairs[0,wh]**2),}

        locbeta=betatriples[:,twh]
        loceta=etatriples[:,twh]
        bbv = np.var(locbeta)
        def Hcov(x, C,):
            djk = np.zeros(x.shape[1])
            djk[x[1] == x[2]] = 1
            return bv * (-C * x[1] * x[2] + djk)

        def qcov(x, C,):
            djk = np.zeros(x.shape[1])
            djk[x[1] == x[2]] = 1
            return bv * (C  + djk)

        p0fit=(0,)
        z= np.product(locbeta[:2],axis=0)
        # zerr=(np.sum(locbeta[::-1]*betaerrtriples[:,twh],axis=0) + np.product(betaerrtriples[:,twh],axis=0))/(np.sum(locbeta[::-1],axis=0) + 1)
        zerr=np.max(betaerrtriples[:2,twh],axis=0)
        popt, pcov = scipy.optimize.curve_fit(Hcov, loceta, z, p0=p0fit,sigma=zerr)
        qopt, qcov = scipy.optimize.curve_fit(qcov, loceta, z, p0=p0fit,sigma=zerr)
        correl=qopt[0]
        # test=linregress(np.product(etatriples[1:],axis=0), -z )
        #

        ntests=np.sqrt(len(set(locbeta[0]) ))
        meanetajk=np.mean(loceta[1]*loceta[2])
        dic['C'] =C= popt[0]
        dic['correltest']=-C*meanetajk
        # dic['C2']=test[0]/bv
        dic['Cstd']=pcov[0,0]**.5
        dic['Cstderr']=pcov[0,0]**.5/ntests
        dic['Crelstd'] = pcov[0, 0] ** .5 / np.abs(popt[0])
        dic['correlrelstd'] = qcov[0, 0] ** .5 / np.abs(popt[0])
        dic['correlstd']=qcov[0,0]**.5

        #Other correlations
        offwh=np.where(loceta[1]!=loceta[2])[0]
        pgopt, pgcov = scipy.optimize.curve_fit(Hcov,loceta[:,offwh], np.product(locbeta[[0,3]],axis=0)[offwh],  p0=p0fit,
                                                sigma=np.max(betaerrtriples[[0,3]][:,twh],axis=0)[offwh] )
        pggopt, pggcov = scipy.optimize.curve_fit(Hcov, loceta[:,offwh],np.product(locbeta[2:],axis=0)[offwh],  p0=p0fit,
                                                  sigma=np.max(betaerrtriples[2:][:,twh],axis=0)[offwh])
        dic['Cgam'],dic['Cgamstderr']=pgopt[0],pgcov[0,0]**.5/ntests
        dic['Cgam2'],dic['Cgam2stderr']=pggopt[0],pggcov[0,0]**.5/ntests
        # code_debugger()
        dic.update({
             'correl':correl,'meanetajk':meanetajk})#, 'tetas':etatriples[:,twh],'tbetas':betatriples[:,twh] })
        # legends.append()
        # slopes.append(popt[1]),consts.append(popt[0]),nonempty.append( np.sqrt(bin1*bin2) ),errslopes.append(pcov[1,1]**.5)
        dic['bv_estim']=np.array((bv,bbv))#,popt[1],qopt[1],pgopt[1],pggopt[1]))
        table.append(dic)
    df=pd.DataFrame(table)
    nonempty,slopes,errslopes,consts=[df[x].values for x in  ['meanetai','slope','errslope','const']]
    if plots=='all':
        plt.legend(legends)
    # plt.show()
    popt=None
    try:
        popt, pcov = scipy.optimize.curve_fit(lambda x,b,a:b+a*x, nonempty, slopes, p0=(0, 0))
        A=-popt[1]
    except:
        print 'ERROR hebbian_slope, bins:',nonempty, slopes
        A=0
    def avgsum(x):
        return np.mean(np.sum(x, axis=0))
    S=eta.shape[0]
    Atheo = 1. / (avgsum(eta ** 2) * (S - 1.) / S)
    def atl(x):
        if len(x.shape)<2:
            return x.reshape(x.shape+(1,))
        return x
    def atsum(x):
        return np.atleast_1d(np.sum(x,axis=0))
    bm=np.mean(consts)  # IMPROVES FIT ON TEST DATA COMPARED TO USING BM FROM CONV (see test_hebbian)

    etamf=1 - bm * (atsum(eta) - atl(df['meanetai'].values))
    Atheos =np.mean( (etamf - atl(df['meanetai'].values)) / (atsum(eta ** 2) - atl(df['meanetai2'].values)) ,axis=1)

    if plots=='A':
        error_on_coef= 1./avgsum(eta)**4 * np.sqrt(avgsum( 4* etavar * eta**2 )  )
        plt.plot(nonempty,Atheos,lw=2)# ,linestyle='--')
        plt.errorbar(nonempty,slopes,yerr=df['stderrslopes'], fmt='o',capsize=5,color='k')
        # plt.plot(nonempty,popt[0]+popt[1]*nonempty,linestyle='--')
        plt.axhline(0,linestyle='--',color='r')

    Ctheos= np.mean(1 / (atsum(eta ** 2) - atl(df['meanetai2'].values)), axis=1)
    correltheos=- Ctheos*df['meanetajk'].values
    if plots=='C':
        plt.plot(nonempty, -Ctheos, lw=2)
        plt.errorbar(nonempty, -df['C'],yerr=df['Cstderr'], fmt='o',capsize=5,color='k')
        plt.axhline(0,linestyle='--',color='r')

    result= {'bins':nonempty, 'Ctheos':Ctheos,'Ctheo': np.mean(Ctheos),'bm_from_fit':np.mean(consts),'bm_from_conv':oldbm,
            'slopes':slopes,'errslopes':errslopes, 'Aslope':A,'Aslope/theo':A/Atheo, 'Atheo':Atheo,'Atheos': Atheos, 'correltheos':correltheos,
            }

    for x in ['meanetajk','correl','C','Cstd','correlstd','Crelstd','correlrelstd','binstart','binend','correltest','stderrslopes','Cstderr',
              'Cgam','Cgam2','Cgamstderr','Cgam2stderr','bv_estim']:
        result[x]=df[x].values
    result['loglike']=hebbian_loglike(beta,bm,bv**.5)
    def dist(a,b):
        # return np.mean(a-b)
        from numpy.linalg import norm
        # return norm(a-b)
        return np.corrcoef(a,b)[0,1]#np.mean( 2*np.abs(a-b)/(a+b) )
    result['Adist']=np.mean((slopes))#dist(slopes,Atheos)
    result['Cdist']=np.mean(result['C'])#dist(result['C'],Ctheos)
    result['betadist']=dist((offdiag(beta)),(offdiag(betatheo)))
    if np.sum(eta.shape) > 1:
        if len(eta.shape)>1:
            result['Arank'] = np.mean([np.corrcoef(getranks(offdiag(beta)), getranks(offdiag(np.multiply.outer(et, et))))[0, 1] for et in list(eta.T)])
        else:
            result['Arank'] = np.corrcoef(getranks(offdiag(beta)), getranks(offdiag(np.multiply.outer(eta, eta))))[0, 1]
    else:
        result['betarank']=0
    return result


def create_synthetic_network_LV_gam0(N,mean_alpha,std_alpha,val=1,return_inv=0):
    import numpy.linalg as la
    S=len(N)
    mat=np.eye(S)
    rows=[]
    for k in range(S):
        notk=range(S)
        notk.remove(k)
        Nmi=N[notk]
        rn=np.random.normal(0,1,(S-2,S-2))
        tmp= np.concatenate([np.atleast_2d(Nmi[1:]),rn])
        AA= np.concatenate([np.atleast_2d(Nmi), tmp.T ])
        QQ= la.qr(AA)[0].T
        QQ=QQ/np.sign(QQ[0,0])
        # Sample x
        x=np.dot(QQ, np.ones(S-1)*mean_alpha)
        x[0]=(val-N[k])/np.sqrt(np.sum(Nmi**2))
        x[1:]+= std_alpha*np.random.normal(0,1,S-2)
        # Rotate
        rows.append(QQ.T)
        row=np.dot(QQ.T,x)
        mat[k,notk]=row
    if val ==1:
        assert np.allclose( N, np.dot( la.inv(mat),np.ones(S) ) )
    if return_inv:
        return mat, rows
    return mat

def invert_QQ(mat,rows,mean_alpha):
    import numpy.linalg as la
    S=len(rows)
    inv=[]
    for k in range(S):
        notk=range(S)
        notk.remove(k)
        inv.append(np.dot(la.inv(rows[k]), mat[k,notk])-np.dot(rows[k].T ,np.ones(S-1)*mean_alpha ) )
    return np.array(inv)[:,1:]


