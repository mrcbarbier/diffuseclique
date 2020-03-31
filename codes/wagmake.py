from wagutils import *
import itertools
from statsmodels.nonparametric.smoothers_lowess import lowess
import pickle
from json import dump,load
import scipy.linalg as la




def reldist_type2(x, y):
    xm, ym = np.mean(x), np.mean(y)
    slope = np.mean((x - xm) * (y - ym) ** 2) / np.mean((x - xm) ** 2 * (y - ym))
    return [slope, 1. / slope][int(np.abs(slope) > np.abs(1. / slope))]


def reldist_odr(x, y):
    import scipy.odr as odr
    def f(B, xx):
        return B[0] * xx + B[1]

    linear = odr.Model(f)
    data = odr.Data(x, y, wd=1, we=1)
    res = odr.ODR(data, linear, beta0=[.4, .0])
    res = res.run()
    b = res.beta[0]
    b = np.array([b, 1 / b])[np.argmin(np.abs([b, 1 / b]))]
    # print reldist_type2(x,y),b
    return b


def reldist(x, y, boot=0, typ='', split=0, strip=1,**kwargs):
    #
    x,y=np.array(x),np.array(y)
    idx = list(np.where((~np.isnan(x))&(~np.isnan(y)))[0])
    if boot:
        idx = list(np.random.choice(idx, replace=True, size=x.size))
    x, y = x[idx], y[idx]
    #
    if strip:
        idx = np.argsort(x)
        z = int(np.floor(len(idx) / 20))
        idx = idx[z:-z]
        x, y = x[idx], y[idx]
        idx = np.argsort(y)
        return reldist(x[idx[z:-z]], y[idx[z:-z]], boot=0, strip=0, typ=typ,**kwargs)
    if split:
        idx = np.argsort(x)
        return np.mean([reldist(x[idx], y[idx], boot=0, split=0, typ=typ,**kwargs) for idx in np.array_split(idx, split)])
    #
    if 'odr' in typ:
        return reldist_odr(x, y)
    if 'type2' in typ or (typ is None and 'type2' in sys.argv):
        return reldist_type2(x, y)
    if 'loglike' in typ:
        if not 'var' in kwargs:
            code_debugger()
        v=kwargs.get('var',1)
        if v is 1:
            print 'No var found'
        if len(v.shape)==2:
            try:
                v=v[np.ix_(idx,idx)]
            except:
                pass
            if v.shape[0]==x.shape[0]:
                # print 'Using covariance matrix'
                return  -np.dot((x-y).ravel(),np.dot(la.inv(v ),(x-y).ravel() ))/x.size
        return -np.mean((x - y) ** 2 / v[idx])
    #
    # Relative distance as correlation rescaled by max variance - like correlation but penalizes != scalings
    cov = np.cov(x, y)
    return cov[0, 1] / np.max([cov[0, 0], cov[1, 1]])


def slopedist(x, y, etai, etaj, boot=0, debug=0, per_species=None, linearize=0,**kwargs):
    if boot:
        idx = list(np.random.choice(range(x.size), replace=True, size=x.size))
        x, y = x[idx], y[idx]
    if not per_species is None:
        i, j, species = per_species
        idx = [np.where(i == z)[0] for z in species]
        slopy = [linregress(etaj[i], y[i])[0] if len(i)>2 else np.nan for i in idx]
        slopx = [linregress(etaj[i], x[i])[0] if len(i)>2 else np.nan for i in idx]
    #
    else:
        idx = np.argsort(etai)
        idx = np.array_split(idx, 4)
        slopy = [linregress(etaj[i], y[i])[0] if len(i)>2 else np.nan for i in idx]
        slopx = [linregress(etaj[i], x[i])[0] if len(i)>2 else np.nan for i in idx]
    loc = np.array([np.median(etai[i]) for i in idx])
    #
    slopx,slopy,loc=[np.array(z) for z in ( slopx,slopy,loc)]

    if linearize:
        good=(~np.isnan(slopx))&(~np.isnan(slopy))
        slopx,slopy,loc=slopx[good],slopy[good],loc[good]
        a = np.argsort(loc)
        slopx, slopy = slopx[a], slopy[a]
        sx, ix = linregress(loc, slopx)[:2]
        sy, iy = linregress(loc, slopy)[:2]
        slopx, slopy = loc * sx + ix, loc * sy + iy
    #
    if 'debug' in sys.argv or debug:
        plt.close('all')
        plot(loc, slopy, hold=1)
        scatter(loc, slopx)
        code_debugger()
    kwargs.setdefault('strip',0)
    if kwargs.get('return_all'):
        return slopx,slopy,reldist(slopy, slopx, **kwargs)
    return reldist(slopy, slopx, **kwargs)

def rowdist(x,y,i,**kwargs):
    kwargs.setdefault('strip',0)
    species=kwargs.get('species',np.unique(i))
    meanx=np.array([np.mean(x[i==s]) if (i==s).any() else 0 for s in species])
    meany=np.array([np.mean(y[i==s]) if (i==s).any() else 0 for s in species])
    return reldist(meanx,meany,**kwargs)

def get_species(df):
    return sorted(set(np.concatenate(df['composition'].values)))

def get_path(exp,return_suffix=0):
    suffix = ''
    if 'bug' in sys.argv:
        suffix += '_debug'
    if 'detrendK' in sys.argv:
        suffix += '_detK'
    elif 'detrendblock' in sys.argv:
        suffix += '_detb'
    elif 'detrend' in sys.argv:
        suffix += '_det'
    if 'cheat' in sys.argv:
        suffix += '_cheat'
    if 'bminfer' in sys.argv:
        suffix += '_bminfer'
    path = Path('data/' + exp + suffix)
    path.mkdir()
    if return_suffix:
        return path,suffix
    return path


def hyperplane_light(df,species,**kwargs):
    df=df.copy()
    from numpy.linalg import lstsq, norm as lanorm, inv as lainv
    from scipy.optimize import least_squares
    S = len(species)
    mat = np.zeros((S, S)) #Returned matrix
    for sidx,s in enumerate(species):
        res=None
        notsidx,others=[list(x) for x in zip(*[(o,oth) for o,oth in enumerate(species) if oth!=s])]
        xs = df[df[s] != 0][species].values
        xnomono=df[(df[s]!=0) & (np.max(df[others],axis=1)!=0 ) ]
        if not xnomono.size:
            print 'hyperplane skipping',s,'only present in monoculture'
            mat[sidx,sidx]=-1
            continue
        xsT = xs.T
        def costfeta(y,weights=None):
            yy = -np.ones(S)
            yy[notsidx] = y
            if weights is None:
                return (np.dot(yy, xsT) +1)
            else:
                return (np.dot(yy, np.sum(weights*xsT,axis=1)/np.sum(weights) ) +1)
        res=least_squares(costfeta,-np.zeros(S-1) )
        row=list(res.x)
        row.insert(sidx, -1)
        mat[sidx]=row
        mat[sidx,np.sum(np.abs(xsT),axis=1)==0]=np.nan
    return mat


def hyperplane(df,species,etamode=0,distances=0,use_var=0,missing=0,**kwargs):
    df=df.copy()
    debug=kwargs.get('debug')
    from numpy.linalg import lstsq, norm as lanorm, inv as lainv
    from scipy.optimize import least_squares
    S = len(species)
    mat = np.zeros((S, S)) #Returned matrix
    compmat=np.zeros((S,S)) #To compare between differnet methods
    errmat=np.zeros((S,S)) #Matrix of stderr on coefficients
    table=[]
    sidx=-1
    res=None
    Kmono=kwargs.get("K",None)
    if Kmono is None:
        Kmono = np.array(
            [df[np.logical_and(np.sum(df[species].values > 0, axis=1) == 1, df[s] > 0)][s].mean() for s in
             species])
        Kmono[np.isnan(Kmono)] = 10 ** -10
    for s in species:
        res,res2=None,None
        sidx+=1
        rsquared=0
        notsidx = [z for z in range(S) if z != sidx]
        xs = df[df[s] != 0][species].values
        xnomono=df[(df[s]!=0) & (np.max(df[[oth for oth in species if oth!=s]],axis=1)!=0 ) ]
        if not xnomono.size:
            print 'hyperplane skipping',s,'only present in monoculture'
            mat[sidx,sidx]=-1
            dic={'species':s,'R2':0,'K':10**-10,'Kvar':0 }
            table.append(dic)
            continue
        if etamode==1:
            print 'basic eta mode'
            xs = xs / np.atleast_1d(Kmono)
            xs[:, np.where(np.isnan(Kmono))[0]] = 0
        xsT = xs.T
        weights=np.ones(xs.shape[0])
        if 'weights' in kwargs:
            weights*=[ kwargs['weights'].get(surv,0) for surv in np.sum(xs>0,axis=1)]
        if distances or debug:
            # print 'USING DISTANCES',distances,debug
            dxs = np.concatenate([xs - x for x in xs]).T  # [1:]-xs[0]
            def costf(y,weights=None):
                yy = -np.ones(S)
                yy[notsidx] = y
                return np.dot(yy, dxs)# -Kmono[sidx]
            try:
                res = least_squares(costf,- np.ones(S - 1))
                if kwargs.get('weights',None) and  not np.allclose(weights,1):
                    print 'Weights+distances not implemented yet'
                    res = least_squares(costf, res.x,kwargs={'weights':weights})
            except Exception as e:
                print 'Failed least_squares',e
                code_debugger()
            ai = list(res.x)
            residuals=res.fun
            ai.insert(sidx, -1)
            mat[sidx] = ai
            Ks=-np.dot(ai,xsT)
            rsquared = 1 - np.sum(residuals ** 2) / np.sum((dxs-np.mean(dxs,axis=1).reshape((S,1)) )**2 )
        if (not distances) or debug:
            def costfeta(y,weights=None):
                # return np.dot(y,xsT[notsidx])-xsT[sidx]+ifelse(etamode=='given',1,Kmono[sidx])
                yy = -np.ones(S)
                yy[notsidx] = y
                if weights is None:
                    return (np.dot(yy, xsT) +1)
                else:
                    return (np.dot(yy, np.sum(weights*xsT,axis=1)/np.sum(weights) ) +1)
            def costfnoeta(y,weights=None):
                # return np.dot(y,xsT[notsidx])-xsT[sidx]+ifelse(etamode=='given',1,Kmono[sidx])
                yy = -np.ones(S)
                yy[notsidx] = y[notsidx]
                if weights is None:
                    return np.dot(yy, xsT) +y[sidx]
                else:
                    raise Exception('NOT READY')
                    return (np.dot(yy, np.sum(weights*xsT,axis=1)/np.sum(weights) ) +1)

            Ks = None
            if etamode:
                try:
                    res2=least_squares(costfeta,-np.ones(S-1) )
                except:
                    code_debugger()
                if kwargs.get('weights',None) and  not np.allclose(weights,1):
                    # code_debugger()
                    res2 = least_squares(costfeta, res2.x,kwargs={'weights':weights})
                comparison=list(res2.x)
                residuals=costfeta(res2.x)
            else:
                x0=-np.ones(S)
                x0[sidx]=1
                res2=least_squares(costfnoeta,x0 )
                if kwargs.get('weights',None) and not np.allclose(weights,1):
                    # code_debugger()
                    res2 = least_squares(costfnoeta, res2.x,kwargs={'weights':weights})
                Ks=res2.x[sidx]
                comparison=list(res2.x[notsidx])
                residuals = costfnoeta(res2.x)
            if  use_var:
                xvarT = df[df[s] != 0][[sp + '_var' for sp in species]].values.T
                xvarT[np.isnan(xvarT)]=0
                try:
                    def costd(yy):
                        dd,vv=xsT,np.clip(xvarT,10**-5,None)
                        tmpx=costfeta(yy)
                        tmpv=( np.dot(yy**2, vv[notsidx]) + vv[sidx]  ) ** .5
                        # code_debugger()
                        # tmpv=( np.sum(vv,axis=0) ) ** .5
                        # tmpv=1
                        # print tmpv
                        return tmpx/tmpv
                    tmpres = list(
                        least_squares(costd, -np.ones(S - 1)).x)
                    tmpres2 = list(
                        least_squares(costd, comparison).x)
                    comparison=tmpres2
                    # print costd(np.array(tmpres)),costd(np.array(tmpres2))
                except:
                    print 'Failure'
                    code_debugger()
            # print 'Final',np.sum(costd(np.array(comparison))**2)
            comparison.insert(sidx, -1)
            if Ks is None:
                Ks=-np.dot(comparison,xsT)
            compmat[sidx]=comparison
            rsquared = 1 - np.sum(residuals ** 2) / np.sum((xsT-np.mean(xsT,axis=1).reshape((S,1)) )**2 )
            if np.isnan(rsquared).any():
                code_debugger()

            # rsquared = 1 - np.sum(residuals ** 2) / np.mean() np.var(xsT)#[sidx])
        # if rsquared<0:
        #     code_debugger()
        if debug:
            code_debugger()
        try:
            def makerr(res):
                from scipy.linalg import svd
                tmp, s, VT = svd(res.jac, full_matrices=False)
                threshold = np.finfo(float).eps * max(res.jac.shape) * s[0]
                s = s[s > threshold]
                VT = VT[:s.size]
                pcov = np.dot(VT.T / s ** 2, VT)
                return np.clip(np.diag(pcov)*2*res.cost/ np.clip(res.jac.shape[0] - res.jac.shape[1],10**-10,None),None,100)

            fres=res
            if fres is None:
                fres=res2
            if fres.jac.shape[1]==S:
                errmat[sidx]=makerr(fres)**.5
            else:
                errmat[sidx,[si for si in range(S) if si != sidx]]=makerr(fres)**.5
        except Exception as e:
            print 'ERROR hyperplane:',e
        dic={'species':s,'R2':rsquared,'K':np.mean(Ks),'Kvar':np.var(Ks) }
        table.append(dic)
    tab=pd.DataFrame(table)
    Ks=np.array([tab.set_index('species')['K'].loc[s] for s in species ])
    if not distances:
        mat=compmat
    np.fill_diagonal(errmat,0)
    #
    # DEAL WITH MISSING PAIRS
    missingpairs=[(i,j) for i in species for j in species if not np.max(np.prod(df[[i,j]].values,axis=1))>0 ]
    for i,j in missingpairs:
        mat[species.index(i),species.index(j)]=np.nan
        mat[species.index(j),species.index(i)]=np.nan
    if missing=='mean':
        mat[np.isnan(mat)]=np.mean(nonan(offdiag(mat))  )
    else:
        mat[np.isnan(mat)] = missing
    if etamode:
        # tab['Kdiff']=tab['K']
        tab['K']=Kmono*tab['K']
        # code_debugger()
        beta=mat
        alpha = mat / np.multiply.outer(Ks, 1. / Ks)
    else:
        alpha=mat
        beta=mat*np.multiply.outer(Ks,1./Ks)
    return alpha,beta,tab,errmat/(np.abs(mat)+10**-10)




def correlcalc(etafull,beta,gamma=0,return_all=1,pad=0,rank=0,**kwargs):
    '''Compute plot of prediction versus theory for means and correlations'''
    def bootstrap(x):
        if not np.sum(x.shape)>0:
            return x
        return np.mean(np.random.choice(x, size=x.size))
    beta=beta.copy()
    S=etafull.shape[0]
    etamean = np.array([bootstrap(etafull[i]) for i in range(S)])
    bm, bv, betatheo = hebbian_getbstats(beta, etamean,**kwargs )  # betaerr=ana[i][j].get('beta_relerr', np.ones(beta.shape)))
    if isinstance(gamma,basestring):
        gamma=np.corrcoef(offdiag(beta),offdiag(beta.T))[0,1]
        # print  '    gamma',gamma
    betatheo = bm + (betatheo - bm) + gamma * (betatheo.T - bm)

    arange=np.arange(S)
    mean_i=np.multiply.outer(np.arange(S),np.ones(S)).astype('int')
    mean_j=mean_i.T
    # code_debugger()
    # bet=beta.copy()
    # bet[bet == 0] = np.nan
    # bet[np.abs(bet)>3.6]=np.nan
    betadiff = beta - betatheo
    diag = np.multiply.outer(np.ones(S), np.eye(S))

    def removeself(mat):
        S = mat.shape[0]
        ss = range(S)
        mat2 = [mat[i][np.ix_([s for s in ss if s != i], [s for s in ss if s != i])] for i in range(S)]
        return np.array(mat2)

    empirical = removeself(np.einsum('ij,ik->ijk', betadiff, betadiff))
    var = np.array([np.nanmean(empirical[i][np.eye(S - 1) != 0]) for i in range(S)]).reshape((-1, 1, 1))
    empirical /= var + 10 ** -15
    empirical -= removeself(diag)
    prediction = removeself(
        - np.multiply.outer(1. / (np.sum(etamean ** 2) - etamean ** 2 + 0.0001), np.multiply.outer(etamean, etamean)))

    def ms(x):
        return np.concatenate([np.nanmean(x, axis=(1)), np.nanmean(x, axis=(2))])

    corr_i=np.multiply.outer(arange,np.ones((S-1,S-1)) ).astype('int')
    corr_j=removeself(np.multiply.outer(np.ones(S),np.multiply.outer(arange,np.ones(S)) ).astype('int'))
    corr_k=removeself(np.multiply.outer(np.ones(S),np.multiply.outer(np.ones(S),arange) ).astype('int'))
    # prediction,empirical=prediction[removeself(diag)==0],empirical[removeself(diag)==0]
    # prediction,empirical=ms(prediction),ms(empirical)  #Makes things significantly worse

    if kwargs.get('remove_zeros',1):
        beta[beta==0]=np.nan

    results=dict( [('mean_theo', offdiag(betatheo)), ('mean_emp', offdiag(beta)), ('corr_theo', prediction.ravel()),
                     ('corr_emp', empirical.ravel()),('mean_etai',etamean[list( offdiag(mean_i))] ),('mean_etaj',etamean[list( offdiag(mean_j))] ),
                   ('corr_etai', etamean[list(corr_i.ravel())] )  ] )


    # code_debugger()
    for z in ('corr_i', 'corr_j', 'corr_k', 'mean_i', 'mean_j'):
        val=locals()[z]
        if 'mean' in z:
            val=offdiag(val)
        else:
            val=val.ravel()
        results[z]=val

    if rank:
        results={i:getranks(results[i]) for i in results}
    results['bm']=bm
    results['bv']=bv

    from scipy.stats import sem, linregress, siegelslopes, theilslopes, ttest_1samp
    try:
        summary={v: linregress(results[v+'_theo'][~np.isnan( results[v+'_emp'])], results[v+'_emp'][~np.isnan( results[v+'_emp'])] )[0] for v in ('mean','corr')}
    except Exception as e:
        print e
        summary={}
    if return_all:
        results.update(summary)
        if pad:
            for k in results:
                if 'mean_' in k:
                    results[k]=np.concatenate([results[k], np.ones(len(results['corr_theo'])-len(results[k]) ) *np.nan ])
                # else:

        return results
    return summary







def infer_bm(eta,meanfield=1,nmat=1,S=None,maxtrials=100,resolution=ifelse('bug' in sys.argv,3,15),use_noise=0, **kwargs):
    from numpy.linalg import lstsq, norm as lanorm, inv as lainv
    from scipy.special import erf
    import time
    Salive=eta.shape[0]
    if S is None:
        S=Salive
    tstart=time.time()
    eta=eta[np.argsort(np.mean(eta,axis=1))]
    mneta=np.mean(np.median(eta,axis=1))
    sdeta=np.std(np.median(eta,axis=1))
    phieta=np.mean( eta>0) *Salive*1./S

    if eta.shape[1] == 1  or Salive<2:
        covmat = np.eye(3)
    else:
        var_mneta = np.array(np.mean([np.random.choice(eta[i], size=maxtrials) for i in range(Salive)], axis=0))
        var_sdeta = np.array(np.std([np.random.choice(eta[i], size=maxtrials) for i in range(Salive)], axis=0))
        var_phieta = np.array(np.mean([np.random.choice((eta[i] > 0), size=maxtrials) for i in range(Salive)], axis=0) *
                              (Salive - np.random.randint(0, 2, maxtrials)) * 1. / S)
        covmat = np.cov([var_mneta, var_sdeta, var_phieta])

    etavec=np.mean(eta,axis=1)
    vare=np.mean(np.var(eta,axis=1)/etavec**1.5)
    # bm_mf = (1. / mneta - 1) / S
    bm_surv=bm=hebbian_stablebm(etavec)
    bs_surv=bs= np.sqrt( (1- np.mean(etavec)**2/np.mean(etavec**2)) /S)
    tab,learntab =None,None
    gamma=0
    if not meanfield:
        #
        if 'table' in kwargs and not kwargs['table'] is None:
            tab = kwargs['table'].copy()
        else:
            def make_mats(bm,bs,gamma):
                etas=[]
                trial=0
                while len(etas)<nmat and trial<maxtrials:
                    trial+=1
                    mat= -genpool(S,mu=bm,sigma=bs,gamma=gamma)[-1]
                    np.fill_diagonal(mat,1)
                    e=np.dot(lainv(mat),np.ones(S))
                    a=np.argsort(e)
                    mat=mat[np.ix_(a,a)]
                    if (e<0).any():
                        e=np.clip(e[a],0,None)
                        e=dosimu(-mat,np.ones(S),np.ones(S),tmax=100,noise=0,x0=e+.001)[-1][-1]
                        e[e<10**-5]=0
                        # print e
                    if vare>0 and use_noise:
                        noise=np.array([np.random.gamma(1./(vare*ee**1.5),(vare*ee**1.5) ) if ee >0 else 0 for ee in e  ])
                    else:
                        noise=1
                    etas.append(e*noise )
                return etas
            #
            learntab=[]
            print 'CREATING TABLE FOR BM INFERENCE'


            if 'nogamma' in sys.argv:
                gammas=[0]
                gammares=1
            else:
                gammas=[0,0.3,0.6]
                gammares=9
            fineres=3*resolution
            for ix,x in enumerate(np.linspace(bm_surv*.7,2*bm_surv,resolution)):
                print '{}/{}'.format(ix+1,resolution)
                for y in  np.linspace(bs_surv*.7,2.4*bs_surv,resolution):
                    for g in gammas:
                        etas=make_mats(x,y,g)
                        mns=[np.mean(e) for e in etas]
                        sds=[np.std(e) for e in etas]
                        phi=[np.mean(e>0) for e in etas]
                        learntab.append({'bm':x,'bs':y,'gamma':g,'mn':np.mean(mns),'sd':np.mean(sds),'phi':np.mean(phi),},)
            learntab=pd.DataFrame(learntab)
        #
            XY=learntab[['bm','bs','gamma']].values
            from sklearn.kernel_ridge import KernelRidge
            from sklearn.model_selection import GridSearchCV
            # ax.plot_wireframe(X, Y, Z, rstride=2, cstride=2,color='k',alpha=.7)
            clf = GridSearchCV(KernelRidge(kernel='rbf', gamma=0.1), cv=5,
                          param_grid={"alpha": [1e0, 0.1, 1e-2, 1e-3],
                                      "gamma": np.logspace(-2, 2, 5)})
            XYsmooth=np.array([x for x in itertools.product( np.linspace(XY[0,0],XY[-1,0],fineres),
                                               np.linspace(XY[0,1],XY[-1,1],fineres) , np.linspace(gammas[0],gammas[-1],gammares) ) ])
            tab = pd.DataFrame({'bm': XYsmooth[:, 0], 'bs': XYsmooth[:, 1], 'gamma': XYsmooth[:, 2]})
            for label in ['mn','sd','phi']:
                Z=learntab[label].values #bandwidth(XY) #kernel='gaussian', bandwidth=.2
                clf.fit(XY,Z)
                Zsmooth = clf.predict(XYsmooth)
                tab[label]=Zsmooth
        for tb in (tab,learntab):
            if tb is None:
                continue
            try:
                dist=[ np.sqrt( np.dot(  (mneta-mn,sdeta-sd,phieta-phi),np.dot(lainv(covmat), (mneta-mn,sdeta-sd,phieta-phi) ) ) )  for mn,sd,phi in tb[['mn','sd','phi']].values]
            except:
                code_debugger()
            # dist = np.abs(mneta*phieta-tb['mn'].values)
            tb['dist']=dist
            tb['p']=[(1-erf( dd/np.min(dist)  )    )/2. for dd in dist]
            # code_debugger()

        if 'bug' in sys.argv:
            for i in ('p','dist','phi','mn'):
                plt.figure()
                if not learntab is None:
                    gammas = learntab['gamma'].unique()
                    plt.subplot(121)
                    plt.imshow(learntab[learntab['gamma']==gammas[0]].pivot('bm', 'bs')[i]), plt.colorbar()
                    plt.subplot(122)
                #
                bms,bss=[['{:.2f}'.format(l) for l in ll]  for ll in (np.sort(tab['bm'].unique()),np.sort(tab['bs'].unique()))]
                gammas = tab['gamma'].unique()
                plt.imshow( tab[tab['gamma']==gammas[0]].pivot('bm', 'bs')[i]), plt.colorbar()
                ax = plt.gca()
                dx,dy=ax.get_xticks(),ax.get_yticks()
                def getidx(lst,z):
                    return [lst[int(i)] if i>=0 and i<len(lst) else '' for i in z]
                ax.set_yticklabels(getidx(bms,dy)),ax.set_xticklabels(getidx(bss,dx)),plt.suptitle(i)
            plt.show()
            # code_debugger()
        bm,bs,gamma=tab[tab['p']>.95 *tab['p'].max()][['bm','bs','gamma']].median()
        tend=time.time()
        # bloc = hebbian_convert(eta, bm, bs, forward=0)
        # bm = bloc['bm_pool']
        print 'Time',tend-tstart, 'bm,bs:',bm_surv,bs_surv, '->',bm,bs,'gamma',gamma#,'bs',bs,'<-',bs_surv
        if np.isnan([bm,bs]).any():
            code_debugger()
    return bm,bs,gamma, tab


def make_groundtruth(S=8,species=None,noise=.08,sizes=(1,2,4,8),replicas=1,nplots=None,plots=None,**kwargs):
    """Create a fake experimental setup to test inference methods."""
    from scipy.misc import comb
    table=[]
    if species is None:
        import string
        ref=string.ascii_lowercase+string.ascii_uppercase+''.join([str(i) for i in range(10)])
        species=[ref[i] for i in range(S)]
    species=np.array(species)
    rs,Ks,beta=genpool(S,**kwargs)
    rs=kwargs.pop('rs',rs)
    Ks=kwargs.pop('Ks',Ks)
    beta=kwargs.pop('beta',beta)
    Aij=beta*np.multiply.outer(rs,1./Ks)
    alpha=beta*np.multiply.outer(Ks,1./Ks)

    def get_true_eq(compo,N=None):
        bet=beta - np.eye(S)
        if isinstance(compo[0],basestring):
            sidx=[list(species).index(i) for i in compo]
        else:
            sidx=compo
        bet=bet[np.ix_(sidx,sidx)]
        eqs = find_eqs(bet, uninvadable=1, stable=1, largest=N is None)
        if not eqs:
            eqs = find_eqs(bet, uninvadable=0, stable=1, largest=N is None)
        if not eqs:
            eqs = find_eqs(bet, uninvadable=0, stable=0, largest=N is None)
        eq= eqs[0]
        if not N is None:
            from scipy.linalg import norm
            eq= eqs[np.argmin([ norm(N-eq)  for eq in eqs   ]) ]
        val=np.zeros(S)
        val[sidx]=eq
        return val
    if plots is None:
        if sizes is None:
            sizes = [2 ** x for x in range(int(np.floor(np.log(S * 1.001) / np.log(2)) + 1))]
            if not S in sizes:
                sizes += [S]
        sizes = np.array(sizes)
        if replicas is None:
            replicas = [int(np.round(S / s)) for s in sizes]
        else:
            replicas = np.atleast_1d(replicas)
            if replicas.shape[0] < sizes.shape[0]:
                replicas = replicas * np.ones(sizes.shape)
        if nplots is None:
            nplots = np.array([min(10, comb(S, s)) if s > 1 else S for s in sizes]) * replicas
        plots=[]
        for size, nrep, npl in zip(sizes,replicas,nplots):
            nsamp=max(1,npl/nrep)
            if npl>comb(S,size):
                samples=list(tuple(x) for x in itertools.combinations(range(int(S)),size))
            elif comb(S,size)<5000:
                allcombs=[tuple(x) for x in itertools.combinations(range(int(S)), size)]
                samples = [allcombs[i] for i in np.random.choice(range(len(allcombs)),int(nsamp),replace=0 )]
            else:
                samples=[tuple(np.random.choice(range(int(S)),size,replace=0)) for i in range(int(nsamp))]
            try:
                nrep = max(1,int(min(nrep, npl / len(samples))))
            except:
                print 'ERROR', size,nrep,npl,nsamp,samples
                code_debugger()
            # print size, nrep,npl,samples, len(samples)==len(set(samples))
            plots+=[species[list(sidx)] for sidx in samples for rep in range(nrep)   ]
    plotn=0
    x0=kwargs.pop('x0',np.ones(len(species)))
    for plot in plots:
        plotn += 1
        sidx=[list(species).index(s) for s in plot]
        print 'Plot {} Species {}'.format(plotn, species[sidx])
        years,results = dosimu(Aij[np.ix_(sidx, sidx)], Ks[sidx], rs[sidx],x0=x0[sidx], noise=noise, evol=1, print_msg=0, **kwargs)
        #print results[-1]/Ks[sidx]
        for year, res in zip(years,results):
            total = np.sum(res)
            dic = {'plot': plotn, 'total': total, 'total/m2': total, 'year': int(np.round(year)), 'richness': len(plot),
                   'composition': tuple(species[sidx]),#'equilibrium':get_true_eq(sidx,plotn) 
                   }
            basekeys=sorted(dic.keys())
            abund = np.zeros(len(species))
            # code_debugger()
            abund[sidx] = res
            dic.update({s: a for s, a in zip(species, abund)})
            table.append(dic)
    df=pd.DataFrame(table)
    df=df[list(basekeys)+list(species)]
    ground={}
    ground.update(kwargs)
    ground.update({'A':Aij-np.diag(rs/Ks),'beta':beta-np.eye(S),'alpha':alpha-np.eye(S),'Ks':Ks, 'rs':rs,'noise':noise,
                   'equilibrium':{compo:get_true_eq(compo) for compo in set(df['composition'].values) }  }
                   )
    return df,ground
