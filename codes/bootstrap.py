from wagmake import *

# OPTIONS IN SYS.ARGV:
#   bug:  DEBUG
# For bootstrap
#   rerun
#   monosplit: split monocultures into two halves, one for infer and one for test
#   noboot: instead of bootstrap, *choice* (does not reduce variance)
#   nogamma: do not search over gamma!=0
#   mf: use meanfield estimate of avg beta
# For show_bootstrap
#   reshow
#   bm: subtract mean from groups of mean_emp and mean_theo  values
#   noperm: do not permute infer and test sets of different bootstrapped replicas (warning: anti-correlations due to high-K monocultures being either in infer or in test)





class Partition():
    def __init__(self,df,criterion,split_mono=1,filter_by={}):
        self.errors=[]
        for k,v in filter_by.get('less',{}).items() :
            df=df[df[k]<v]
        for k, v in filter_by.get('more', {}).items():
            df = df[df[k] > v]
        for k,v in filter_by.get('equal',{}).items() :
            df=df[df[k]==v]

        mono=df[df['richness']==1]
        monoplots=list(mono['plot'].unique() )
        poly=df[(df['richness']>1)&(~df['plot'].isin(monoplots) )]
        if criterion:
            op, key, val = criterion
            nomono=not mono[op(mono[key],val)].shape[0]
        else:
            nomono=False

        infer,test,monotest=[],[],[]
        if split_mono and not nomono :
            for s,gp in mono.groupby('composition'):
                repeats=np.unique(gp['plot'])
                if len(repeats)>1:
                    repeats=repeats[np.random.permutation(range(len(repeats)))]
                    mid=np.ceil(len(repeats)/2.).astype('int')
                    infer+=list(repeats[:mid])
                    monotest+=list(repeats[mid:])
                else:
                    self.errors.append('No mono split for {}'.format(s)  )
                    infer.append(repeats[0])
                    monotest.append(repeats[0])
        else:
            if not nomono:
                infer+=[pl for pl,gp in mono.groupby('plot')]
            monotest+=[pl for pl,gp in mono.groupby('plot')]

        if not criterion:
            all=[pl for pl,gp in poly.groupby('plot')]
            infer+=all
            test+=[pl for pl,gp in poly[poly['richness']==poly['richness'].max()].groupby('plot')]
        else:
            infer+=[pl for pl,gp in poly[op(poly[key],val)].groupby('plot')]
            test+=[pl for pl,gp in poly[~op(poly[key],val)].groupby('plot')]
        self.infer=df[df['plot'].isin(infer)]
        self.test=df[df['plot'].isin(test)]
        self.monotest=df[df['plot'].isin(monotest)]
        overlaps=set(self.infer['plot']).intersection(self.test['plot']),set(self.infer['plot']).intersection(self.monotest['plot'])
        if overlaps[0] or overlaps[1]:
            print 'Overlap',len(overlaps[0]),len(overlaps[1])
            # code_debugger()
        if self.errors:
            print self.errors


def running_average(dfm,dfc,species,**kwargs):
    from bisect import bisect_left
    nbins = kwargs.get('runbins', 10)
    avgd = []
    for dfb in (dfm, dfc):
        k = [key for key in dfb if 'emp' in key and not 'runavg' in key]
        if len(k) > 1:
            raise Exception('show_bootstrap runavg: multikey not implemented')
        k = k[0]
        xlab = k.replace('emp', 'theo')
        dff=dfb.copy()
        dff['oldindex']=dfb.index
        dff = dff.sort_values(xlab)
        xs = dff[xlab].values
        ys = dff[k].values


        if kwargs.get('old',0):
            """OBSOLETE"""
            subbin = float(kwargs.get('runsub', 3))
            binrg = np.linspace(np.min(xs), np.max(xs) + 10 ** -5, nbins + 1)
            binpos = np.array([bisect_left(list(xs) + [np.inf], b) for b in binrg])
            binlist = zip(binpos[:-1], binpos[1:], np.maximum(1, (binpos[1:] - binpos[:-1]) / subbin).astype('int'))
            binlist = [b for b in binlist if b[0] != b[1]]
            binlist[-1] = tuple(binlist[-1] + np.array((0, 1, 0), dtype='int'))
            binlist = [tuple(b) + (b[2] + min(0, b[0] - b[2]),) for b in binlist]  # w  = binavg+min(0,binstart-binavg)
            yavg = np.concatenate(
                [np.convolve(ys[max(0, binstart - binavg):min(len(xs), binend + binavg)], np.ones(binavg) * 1. / binavg,
                             'same')[w:w + (binend - binstart)] for (binstart, binend, binavg, w) in binlist])
        else:
            binwid=(kwargs.get('binmax',np.percentile(xs,90))-kwargs.get('binmin',np.percentile(xs,10)))/nbins
            if k+'_runavg' in dff:
                print 'processing...'
            left=np.searchsorted(xs,xs-binwid)
            right=np.searchsorted(xs,xs+binwid)
            yavg = np.array([np.nanmean(ys[x1:x2])  for x1,x2 in zip(left,right)   ])
            if k+'_runavg' in dff:
                nrep=len(dff['replica'].unique())*1.
                # Standard error wiht average number of points *per bootstrap replica*
                ystd= np.array([ np.std(nonan(ys[x1:x2]))/np.sqrt((x2-x1)/nrep)  for x1,x2 in zip(left,right)   ])
                dff[k+'_05'],dff[k+'_95']=-1*ystd,1*ystd
            if k+'_runavg' in dff:
                print 'processed'

        dff[k] = yavg
        dff=dff.sort_values('oldindex')
        dff.drop('oldindex',axis=1)
        avgd.append(dff)

    adfm, adfc = avgd
    def indices(c):
        i = np.multiply.outer(range(len(c)), np.ones(len(c))).astype('int')
        return zip(offdiag(i), offdiag(i.T))


    """ RECONSTRUCT MATRICES """
    agg=[np.mean,np.median][1]

    sadfm = adfm.sort_index()
    sadfm['pair']=[(x,y) for x,y in sadfm[['mean_i','mean_j']].values]
    present = set(sadfm['pair'].values)
    runemp = [agg(sadfm.loc[sadfm['pair'] == pair,'mean_emp'].values) if pair in present else 0 for pair in
              (indices(species))]
    runtheo = [agg(sadfm.loc[sadfm['pair'] == pair,'mean_theo'].values) if pair in present else 0 for pair in
               (indices(species))]

    rndfm=sadfm[['pair','mean_emp']].copy()
    rndfm['pair']=rndfm['pair'].values[list(np.random.permutation(range(rndfm.shape[0])) )]
    runrnd = [agg(rndfm.loc[rndfm['pair'] == pair,'mean_emp'].values) if pair in present else 0 for pair in
               (indices(species))]
    # testheo=np.mean(results['mean_theo'].values,axis=0)
    # assert np.allclose(runtheo,testheo)

    S=len(species)
    runavgmat = np.zeros((S, S))
    runavgmat[np.eye(S) == 0] = runemp
    theomat = np.zeros((S, S))
    theomat[np.eye(S) == 0] = runtheo
    rndmat = np.zeros((S, S))
    rndmat[np.eye(S) == 0] = runrnd
    # code_debugger()

    return adfm,adfc,runavgmat,theomat,rndmat




def cc_analysis(bootdata,species,N,Kmono,eta,beta=None,betatheo=None,bm=None,all_plots=1,betavar=None,runavgmat=None,**kwargs):
    incomp=np.array([[ (s in c)  for s in species] for c in eta.composition])
    # incomp=1.

    figs={}
    measures={}
    def jaillard(eta,beta,K,bm=None):
        if bm is None:
            bm=np.mean(offdiag(beta))
        compsize=[np.array([len(c)*1. for c in eta.composition if s in c]) for s in species]
        invs=[np.mean(1/S) for S in compsize]
        compo_effect=1+(K/np.mean(K)-1)*invs
        bb=beta.copy()
        np.fill_diagonal(bb,np.nan)
        bi=np.array([np.nanmean(bb[:,i]) for i in range(bb.shape[0])])
        kr,br=K/np.mean(K), (1-bi)/(1-bm)
        ie= (kr-1)*(1-br)*[np.mean((S-1)/(k+S-1)/(1+(S-1)*b ) ) for S,k,b in zip(compsize,kr,br)]
        inter_effect=1./bm*(ie+1)
        return  compo_effect, inter_effect

    if beta is None:
        beta = -hyperplane_light(eta, species)

    if all_plots:
        xs, ys=jaillard(eta,beta,Kmono)
        figs['jaillard']=plt.figure()
        scatter(xs,ys,hold=1)
        for s,x,y in zip(species,xs,ys):
            plt.annotate(s,(x,y))

    beta_compar=None
    sumratio=(np.sum(eta[species].values,axis=1)/np.sum(eta[species].values*incomp,axis=1)  ).reshape((-1,1))
    eta2=eta.copy()
    eta2[species]=(eta[species].values*incomp)
    eta2[species].values[eta2[species]<0.01]=0
    beta2 = -hyperplane_light(eta2, species)
    #beta2[np.isnan(beta2)]=0

    alpha, beta_from_alpha, alpha_details, alpha_relerr=hyperplane(N, species, etamode=0)
    obeta, obeta_relerr=hyperplane(eta, species, etamode='given')[0::3]
    alpha,obeta=-alpha,-obeta
    Khyp=alpha_details["K"]

    if all_plots:
        figs['clean_comparison'] =plt.figure()
        plt.hist(nonan(np.concatenate([np.dot(beta, e) for e in eta[species].values])), bins=50, histtype='step',density=1)
        plt.hist(nonan(np.concatenate([np.dot(beta2, e) for e in eta2[species].values])), bins=100, histtype='step',density=1)
        plt.hist(nonan(np.concatenate([np.dot(alpha, n)/Khyp for n in N[species].values])), bins=100, histtype='step',density=1)
        plt.legend(['beta.eta','beta_clean.eta_clean','alpha.N / K'])

    octoeta = eta[[len(c) >= 4 for c in eta.composition]][species].mean()

    tests=[]
    trueoctoeta,truebeta=octoeta,beta
    def mkslope(mat,eta,rem0=kwargs.get('remove_zeros',0)):
        eta,mat=np.array(eta),np.array(mat)
        res=[]
        for i in range(len(eta)):
            noti=range(len(eta))
            noti.remove(i)
            xs=eta[noti]
            ys=mat[i,noti]
            xs,ys=np.array(xs),np.array(ys)
            good = (~np.isnan(xs))&(~np.isnan(ys))
            xs, ys = xs[good], ys[good]
            if rem0:
                good=(ys!=0)
                xs,ys=xs[good],ys[good]
            if len(xs)>2:
                slope, intercept, r, p, stderr = linregress(xs, ys)
            else:
                slope=np.nan
            # print i,xs, ys ,slope
            res.append(slope)
        return np.array(res)

    testidx=0
    for octoeta,beta in [(trueoctoeta,truebeta),]+tests:
        testidx+=1
        if 'betalim' in sys.argv:
            beta=np.clip(beta,-2,2)
        if kwargs.get('remove_err', 1) and 'bs' in kwargs:
            beta[betavar>kwargs['bs']**2]=np.nan
        octosort=np.argsort(octoeta).values
        figs['shapecomp_{}'.format(testidx)] =plt.figure(figsize=np.array(mpfig.rcParams['figure.figsize'])*(1.3,2.5))
        if not 'bminfer' in sys.argv :
            tit=kwargs.get('title',r'')
            plt.suptitle(tit.capitalize() )#+ r' ($\bar \beta$ obtained from $\beta$)')
        else:
            plt.suptitle(kwargs.get('title',r'')+r' ($\bar \beta$ obtained from $\eta$)')
        if betatheo is None:
            bm, bv, betatheo = hebbian_getbstats(beta,octoeta.values,get_bm_from_beta= not 'bminfer' in sys.argv )

        def good(mat):
            mat0=offdiag(mat)
            return mat0[~np.isnan(mat0)]

        if not 'bs' in kwargs:
            bs=(np.var(good(beta-betatheo)) - np.mean(good(betavar) ))**.5
            if np.isnan(bs):
                print '\n  BS ERROR!',np.var(good(beta-betatheo)),np.mean(good(betavar)),'\n'
                bs=np.var(good(snap_beta(beta,octoeta)-betatheo ))**.5
            bs2=np.var(good(snap_beta(beta,octoeta)-betatheo ))**.5
            # code_debugger()
            print '   VARIANCES', np.var(good(beta)) ** .5, np.var(good(beta - betatheo)) ** .5, bs2, bs
        else:
            bs=kwargs['bs']
        compmats=[create_synthetic_network_LV_gam0(octoeta, bm, bs) for n in range(ifelse('all' in sys.argv,5000,500))]
        if kwargs.get('use_external_bm',0):
            empbm = bm  # np.mean(offdiag(beta))
            betatheo=np.median([b for  b in compmats],axis=0)
        else:
            #todo
            empbm=bm

        if not good(beta).any():
            continue
        vmax=np.max(np.abs(good(beta)-empbm))
        plt.subplot(321),plt.title(r'Empirical $\beta_{ij}$')
        plt.colorbar(plt.imshow(beta[np.ix_(octosort, octosort)],vmin=empbm-vmax,vmax=bm+vmax,cmap='seismic'))
        plt.yticks(range(len(species)),[s[:3] for s in np.array(species)[octosort]], )
        plt.xticks([])
        # plt.gca().tick_params(labeltop=1,labelbottom=0)
        #
        plt.subplot(322),plt.title(r'Theoretical $E[\beta_{ij}|\eta_i,\eta_j]$')
        # vmax=np.max(np.abs(good(betatheo)-empbm))
        betatheoshow=betatheo.copy()
        betatheoshow[betatheo==0]=bm
        plt.colorbar(plt.imshow(betatheoshow[np.ix_(octosort, octosort)],vmin=empbm-vmax,vmax=empbm+vmax,cmap='seismic'))
        plt.yticks([])#range(len(species)),[s[:3] for s in np.array(species)[octosort]], )
        plt.xticks([])
        etasort = np.array([octoeta[s] for s in octosort])

        alphatheo=betatheo*np.multiply.outer(Kmono,1./Kmono)

        def offmean(mat,row=0):
            tmp=mat.copy()
            np.fill_diagonal(tmp,np.nan)
            if row:
                return [np.nanmean(tmp[s, :]) for s in octosort]
            return np.array([np.nanmean(tmp[:,s]) for s in octosort ])

        titles=['column mean','column slope','row mean','row slope',]
        theo=[ offmean(betatheo), mkslope(betatheo.T,octoeta )[octosort], # offmean(alphatheo),
               offmean(betatheo,row=1),  mkslope(betatheo,octoeta )[octosort],#offmean(alphatheo,row=1),
         ]
        def interval(x,sd=False):
            # return np.mean(
            if sd:
                return    np.array((np.median(x,axis=0)-np.percentile(x,16,axis=0),np.percentile(x,84,axis=0)-np.median(x,axis=0) ))
                # return    np.array((np.median(x,axis=0)-np.percentile(x,5,axis=0),np.percentile(x,95,axis=0)-np.median(x,axis=0) ))
            return    np.array((np.median(x,axis=0)-np.percentile(x,5,axis=0),np.percentile(x,95,axis=0)-np.median(x,axis=0) ))
        width=[ interval([offmean(b) for b in compmats  ] ) , interval([mkslope(b.T,octoeta )[octosort] for b in compmats  ] ),# None,
                interval([offmean(b,row=1) for b in compmats]),  interval([mkslope(b,octoeta )[octosort] for b in compmats  ] ), #None,

                ]
        widthsd = [interval([offmean(b) for b in compmats],sd=1),
                 interval([mkslope(b.T, octoeta)[octosort] for b in compmats],sd=1),# None,
                 interval([offmean(b, row=1) for b in compmats],sd=1),
                 interval([mkslope(b, octoeta)[octosort] for b in compmats],sd=1), #None,

                 ]

        if not runavgmat is None:
            runavg=[None,mkslope(runavgmat.T,octoeta)[octosort],None,mkslope(runavgmat,octoeta)[octosort],]
        else:
            runavg=[None,None,None,None,]

        slope,slopeT=mkslope(beta.T,octoeta)[octosort],mkslope(beta,octoeta)[octosort]
        if np.isnan(slope).any() or np.isnan(slopeT).any():
            print '    NaN in slope'
            slope[np.isnan(slope)]=0
            slopeT[np.isnan(slopeT)]=0
            # code_debugger()
        colmean,rowmean=offmean(beta),offmean(beta.T)

        reldists={0:[reldist(offmean(b ),theo[0] ,strip=0) for b in compmats  ] ,
                  1:[reldist(mkslope(b.T,octoeta )[octosort],theo[1] ,strip=0) for b in compmats  ] ,
                  2:[reldist(offmean(b.T ),theo[2] ,strip=0) for b in compmats  ] ,
                  3:[reldist(mkslope(b,octoeta )[octosort],theo[3] ,strip=0) for b in compmats  ]
                  }

        for idm,mat in enumerate([colmean,slope,rowmean,slopeT]):
            #
            if mat is None:
                continue
            if 'Alpha' in titles[idm]:
                continue
            plt.subplot(323+idm)

            xs,thys=np.array(etasort),theo[idm]
            ys=thys
            slope, intercept, r, p, stderr = linregress(xs, ys)
            wid=width[idm]
            widsd=widthsd[idm]
            from datatools import errorfill
            exs=np.median(xs)+1.3*(xs-np.median(xs))
            # errorfill(exs,intercept+slope*exs,yerr=wid,hold=1,)
            plt.plot(xs,ys,color='k',lw=2)
            plt.plot(xs,ys+wid[1],color='k',linestyle='--')
            plt.plot(xs, ys-wid[0],color='k',linestyle='--')
            plt.fill_between(xs,ys - widsd[0], ys + widsd[1], color='g',alpha=.3)
            # plt.plot(xs,  color='k', linestyle='--')

            # plt.colorbar(plt.imshow(tmp,vmin=bm-vmax,vmax=bm+vmax,cmap='seismic'))
            xs,ys=np.array([octoeta[s] for s in octosort]),mat

            if kwargs.get('detailed_labels',0) or 'debug' in sys.argv:
                EvsA=""
                if not runavg[idm] is None:
                    EvsA=', vs A {:.2g}'.format(reldist(thys, runavg[idm],strip=0))
                plt.title(titles[idm].capitalize() + ' E vs beta {:.2g}{}'.format(reldist(thys,ys,strip=0),EvsA ) )
                data=pd.DataFrame({'x':xs,'y':ys} )
                sns.regplot('x','y',data,line_kws={'alpha':.1},ci=68)
            else:
                plt.scatter(xs,ys,c='k')
                plt.title(titles[idm].capitalize())

            for s, x, y in zip(np.array(species)[octosort], xs, ys):
                plt.annotate(s[:3], (x, y))

            # plt.scatter(xs, ys, c='k')
            slope, intercept, r, p, stderr = linregress(xs, ys)
            # plt.plot(xs,intercept+slope*xs,color='k')
            import seaborn as sns
            plt.gca().set_aspect('auto')
            plt.xlim(xmin=exs[0],xmax=exs[-1])
            plt.xlabel(''),plt.ylabel('')

            rdic={'eta':xs,'slope':ys,'theo':theo[idm],'reldist':reldist(ys,theo[idm],strip=0) }
            for perc in (5,50,95):
                rdic['reldist_gen{}'.format(perc)]=np.percentile(reldists[idm],perc)
            rdic['var']=(wid[1]+wid[0]) **2
            rdic['loglike']=-np.mean( (  (theo[idm]-ys  )**2 /rdic['var'] ) )

            sij = np.arange(len(octosort)).astype('int')
            i = np.multiply.outer(sij, np.ones(len(octosort)).astype('int') )
            i, j = offdiag(i), offdiag(i.T)
            xs=np.array(xs)

            print '     LOG LIKELIHOOD',titles[idm],rdic['loglike']
            # code_debugger()
            if np.isnan(rdic['reldist']):
                code_debugger()


            remnan=(~np.isnan(offdiag(beta[np.ix_(octosort,octosort)])))
            i,j=i[remnan],j[remnan]

            def sort(m):
                return offdiag(m[np.ix_(octosort,octosort)])

            def correct(m):
                return sort(m)[remnan]

            if 'mean' in titles[idm]:
                if 'col' in titles[idm]:
                    measures['colmean'] = rdic
                if 'row' in titles[idm]:
                    measures['rowmean'] = rdic

            if 'slope' in titles[idm]:
                if 'col' in titles[idm]:
                    # print ys,theo[idm]
                    slopx, slopy, logl = slopedist(correct(beta), correct(betatheo), xs[[z for z in j ]],
                                                   xs[[z for z in i]], boot=0,
                                                   per_species=(j, i, sij), return_all=1, var=rdic['var'],
                                                   typ='loglike')
                    measures['colslope'] = rdic
                if 'row' in titles[idm]:
                    slopx, slopy, logl = slopedist(correct(beta), correct(betatheo), xs[[z for z in i]],
                                                   xs[[z for z in j]], boot=0,
                                                   per_species=(i, j, sij), return_all=1, var=rdic['var'],
                                                   typ='loglike')
                    measures['rowslope'] = rdic
                    # print sorted(correct(beta)),'\n', sorted(correct(betatheo)), '\n',xs[[z for z in i]],'\n',i
                    # print '       overall',measures['overall']['loglike'],-np.mean(offdiag(beta - betatheo) ** 2 )
                print '         check (fct used in pattern_score):', logl
                # code_debugger()


                if 'col' in titles[idm]:
                    # varb=np.cov(np.array([ correct(b-betatheo) for b in compmats]).T,ddof=0)
                    # varb TURNS OUT NOT TO BE INVERTIBLE, SO I USE THE THEORETICAL PATTERN INSTEAD
                    SS=octoeta.shape[0]
                    etai=np.multiply.outer(octoeta,np.ones(SS ))
                    ii=np.multiply.outer(np.arange(SS),np.ones(SS )).astype('int')
                    varb=np.eye(correct(betatheo).shape[0]) - np.multiply.outer(correct(etai.T),correct(etai.T))/correct(np.sum(etasort**2)-etai**2)
                    varb=varb*bs**2
                    varb2=np.eye(sort(betatheo).shape[0]) - np.multiply.outer(sort(etai.T),sort(etai.T))/sort(np.sum(etasort**2)-etai**2)
                    varb2=varb2*bs**2
                    measures['overall'] = {'loglike': reldist(correct(beta),correct(betatheo), strip=0,typ='loglike',  var=varb),
                                           'loglikesimu':  np.median([reldist(sort(b),sort(betatheo), strip=0,typ='loglike',  var=varb2) for b in compmats]),
                                            'var': varb2,

                                           }
                    print '      loglike {} loglike simu {}'.format(measures['overall']['loglike'],measures['overall']['loglikesimu'])
                    print '      (check uncorr: {})'.format( - np.mean( (correct(beta)-correct(betatheo) )**2/np.diag(varb) ) ), - np.mean( correct(np.tril(betatheo)) **2/np.diag(varb) ), - np.mean( (correct(beta)-correct(betatheo) )**2/bs**2)


    octoeta,beta=trueoctoeta,truebeta
    return figs,measures




def bootstrap(df,filter_by='drop2',scaling='beta',hyperoptions='',split_criterion='<=4',split_mono=1,gamma=None,ntrials=40,
              fill_missing_mono=0, testmean=['all','plot','none'][1],infermean=['all','plot','none'][1],
              **kwargs):
    #beta=hyperplane(usedf, species, etamode=0, distances='dist' in hyperoptions, weights=weights)))
    def bootstrap(x):
        if 'noboot' in sys.argv:
            return np.random.choice(x)
        val=np.random.choice(x, size=x.size)
        if 'logboot' in sys.argv:
            if np.max(val)<=0:
                return 0
            val= 10**np.mean(np.log10(val[val>0]))
            return val
        elif 'quadboot' in sys.argv:
            return np.mean(val**2)**.5
        elif 'medianboot' in sys.argv:
            return np.median(val)
        else:
            return np.mean(val)
        # return np.array([np.mean(np.random.choice(x, size=x.size)) for z in range(x.size)])
    from mpltools import special as pltspecial, style as pltstyle
    species=get_species(df)
    years = sorted( np.unique(df['year'] ))
    df['composition']=[tuple(c) for c in df['composition']]

    # REPLACE MISSING MONOCULTURES BY K hyperplane
    missing_mono= [  s for s in species if not (s,) in set(df['composition'].values ) ]
    if fill_missing_mono and missing_mono:
        print 'GETTING K FROM HYPERPLANE FOR',missing_mono
        Kmono = np.array([df[df['composition'] == (s,)][s].mean() for s in species])
        Kmono[np.isnan(Kmono)] = 0
        for trial in range(ntrials):
            locdf=pd.DataFrame([dict([(k, bootstrap(gp[k])) for k in species]
                                   + [(k, gp[k].values[0]) for k in gp if not k in species] + [('composition', compo)])
                              for compo, gp in df.groupby('composition') for x in range(len(gp))])
            Khyper=hyperplane(locdf, species, etamode=0, distances='dist' in hyperoptions,
                           use_var='var' in hyperoptions, weights={})[2]['K'].values
            for s in missing_mono:
                dic={ss:0 for ss in species}
                dic.update({ 'composition':(s,),'richness':1 , 'plot':'filler_{}'.format(s),'year':0 ,s:np.maximum(.01,Khyper[species.index(s)]) })
                df=pd.concat([df,pd.DataFrame([dic ]  ) ])
    elif scaling=='beta' and missing_mono:
        print 'REMOVING SPECIES WITHOUT MONOCULTURE',missing_mono
        species=[s for s in species if not s in missing_mono]
        df['composition']=[tuple(s for s in c if s in species ) for c in df['composition'].values ]

    # SETTING UP FILTERS FOR PARTITION
    if isinstance(filter_by,basestring):
        if 'drop' in filter_by:
            filter_by={'more':{ 'year':years[int(filter_by.replace('drop',''))-1 ]} }
        else:
            filter_by={}

    Sall=len(species)
    Smax=df['richness'].max()
    if not split_criterion:
        split_criterion={}

    if isinstance(split_criterion,basestring):
        txtcriterion=split_criterion
        split_criterion={'nofull':(np.less,'richness',Smax),'12':(np.less,'richness',3),'<=4':(np.less,'richness',5),'<=8':(np.less,'richness',9),'<=9':(np.less,'richness',10),
        '<=half':(np.less,'richness',Smax/2+1), 'nomono':(np.greater,'richness',1),'>2':(np.greater,'richness',2),'':None,
                         'no16':(np.not_equal,'richness',16)
                         }[split_criterion]
    else:
        txtcriterion=''
    bootdata=[]
    partitions=[]
    for trial in range(ntrials):
        print trial
        partition=Partition(df,split_criterion,split_mono=split_mono and scaling=='beta',filter_by=filter_by)

        # INFERENCE
        weights = {}

        usedf = partition.infer.copy()
        monocultures=[ s for s in species if (s,) in [tuple(x) for x in usedf['composition']]  ]
        Kmono=np.array([bootstrap(usedf[usedf['composition']==(s,)][s]) if s in monocultures else np.nan for s in species ])
        Kmono[np.isnan(Kmono)]=0#np.mean(Kmono[~np.isnan(Kmono)])

        if infermean=='none':
            N=usedf[list(species)+['composition']]
        elif infermean=='plot':
            N=pd.DataFrame( [dict([(k,bootstrap(gp[k]) ) for k in species ]
                                + [(k,gp[k].values[0]) for k in gp if not k in species]+ [('composition',gp['composition'].values[0] )]    )
                           for pl, gp in usedf.groupby('plot') for x in range(len(gp) )   ] )
        elif infermean=='all':
            N=pd.DataFrame( [dict([(k,bootstrap(gp[k]) ) for k in species ]
                                + [(k,gp[k].values[0]) for k in gp if not k in species]+ [('composition',compo)]    )
                           for compo, gp in usedf.groupby('composition') for x in range(len(gp) )   ] )

        eta = N.copy()
        if scaling=='alpha':
            hyp=hyperplane(eta, species, etamode=0, distances='dist' in hyperoptions,use_var='var' in hyperoptions,weights=weights)
            beta=-hyp[0]
            tab=hyp[2]
            Khyper=tab['K'].values
        else:
            eta[species]=eta[species]/(Kmono+10**-15)
            if kwargs.get('light',1):
                beta=-hyperplane_light(eta,species)
                Khyper=np.zeros(len(species))
            else:
                hyp=hyperplane(eta, species,distances='dist' in hyperoptions, use_var='var' in hyperoptions, etamode='given',
                           K=Kmono, weights=weights)
                beta=-hyp[1]
                tab=hyp[2]
                Khyper=tab['K'].values
        if 'betalim' in kwargs:
            beta=np.clip(beta,-kwargs['betalim'],kwargs['betalim'])
            # beta[np.abs(beta) > kwargs['betalim']] = np.nan
        bootdata.append({'Kmono':Kmono,'Khyper':Khyper, 'eta':eta, 'beta':beta,'species':species})
        partitions.append(partition)
    bootdata=pd.DataFrame(bootdata)

    if scaling=='alpha':
        if txtcriterion=='nomono':
            print hyperplane(df[df['richness'] > 1], species, etamode=0, distances='dist' in hyperoptions,
                            use_var='var' in hyperoptions, weights=weights)[2]['K'].values
            monocultures=[ s for s in species if (s,) in [tuple(x) for x in df['composition']]  ]
            Kmono = np.array(
                [bootstrap(df[df['composition'] == (s,)][s]) if s in monocultures else np.nan for s in species])
            bootdata['Kmono']=[Kmono for x in range(bootdata.shape[0])]
            print np.mean(bootdata['Kmono'],axis=0),np.mean(bootdata['Khyper'],axis=0)
            scatter(np.concatenate(bootdata['Kmono']),np.concatenate(bootdata['Khyper']) )

    # GETTING ALL THE ETA's
    etalocs={}
    for trial in range(ntrials):
        usedf = partitions[trial].monotest.copy()
        Kmono=np.array([bootstrap(usedf[usedf['composition']==(s,)][s]) if s in monocultures else np.nan for s in species ])
        Kmono[np.isnan(Kmono)]=0#np.mean(Kmono[~np.isnan(Kmono)])

        truecomps=[]
        for comp,gp in partition.test.groupby('composition'):
            compo=tuple([c for c in comp if gp[c].mean() > 0 and Kmono[species.index(c)] > 0])
            truecomps+=[compo for z in range(gp.shape[0])]
        partition.test['composition']=truecomps

        for compo,gp in partition.test.groupby('composition'):
            if not compo:
                print 'No species in partition',comp
                continue
            sidx = [species.index(c) for c in compo]
            if testmean=='none':
              Nloc=gp[list(compo)].values.T
            elif testmean=='plot':
                Nloc = np.array([[bootstrap(pgp[k]) for k in compo] for z,pgp in gp.groupby('plot') for x in range(len(pgp)) ]).T
            elif testmean=='all':
                Nloc = np.array([[bootstrap(gp[k]) for k in compo]  for x in range(len(gp)) ]).T

            etaloc=Nloc.copy()
            if scaling=='beta':
                etaloc/=Kmono[sidx].reshape((-1,1))
            etalocs.setdefault(compo,{})
            etalocs[compo][trial]=(etaloc)

    Smax=np.max([len(c) for c in df['composition']])
    tabs = {}

    if len(etalocs.values())>5:
        print 'MANY ETALOCS!',len(etalocs.values())
        # code_debugger()
        etalocs={c:etalocs[c] for c in etalocs if len(c)>max([len(z) for z in etalocs])*.5 }
        print 'NOW',len(etalocs.values())

    allbeta = np.nanmean(list(bootdata['beta'].values), axis=0)
    betaerr =np.nanstd(list(bootdata['beta'].values), axis=0)
    betaerr[np.isnan(betaerr)] = np.max(betaerr[~np.isnan(betaerr)])
    betaerr[betaerr == 0] = np.max(betaerr)
    # code_debugger()
    joineta=np.array([ np.random.choice(np.concatenate([ etalocs[c][z][c.index(s)]  for c in etalocs if s in c for z in etalocs[c]]),size=100)
                       if True in [s in c for c in etalocs] else np.zeros(100)
                       for s in species ])
    for compo in etalocs:
        alleta=np.concatenate(etalocs[compo].values(),axis=1)
        if not 'bminfer' in sys.argv:
            bmall, bvall = hebbian_getbstats(allbeta,joineta, get_bm_from_beta=1,betaerr=betaerr)[:2]
            if np.isnan(bmall):
                code_debugger()
            bsall,gamall,tab=bvall**.5,0,None
        else:
            bmall, bsall, gamall, tab = infer_bm(alleta, S=Smax, table=None,use_noise=1,
                                         meanfield=kwargs.get('use_meanfield_bm', 'mf' in sys.argv))

        tabs[compo]=tab
        print  'bm,bs,gamma:',bmall,bsall,gamall, 'mu,sigma',bmall*Smax,bsall*np.sqrt(Smax)

    results=[]
    for trial in range(ntrials):
        print trial
        # TESTING
        usedf = partitions[trial].monotest.copy()
        Kmono=np.array([bootstrap(usedf[usedf['composition']==(s,)][s]) if s in monocultures else np.nan for s in species ])
        Kmono[np.isnan(Kmono)]=0
        beta = bootdata.loc[trial]['beta']
        if kwargs.get('take_mean_beta', 0):
            beta = np.mean(bootdata['beta'], axis=0)

        for comp, gp in partition.test.groupby('composition'):
            compo=tuple([c for c in comp if gp[c].mean()>0 and Kmono[species.index(c)]>0])#0.001*np.median(gp[c].values) ]
            if not compo or not compo in etalocs:
                # print "NO COMPO",compo, comp
                # if compo == comp:
                    # code_debugger()
                continue
            if not trial in etalocs[compo]:
                # print "NO TRIAL",trial
                continue
            # print "======================= OK!", compo
            try:
                sidx = [species.index(c) for c in compo]
                etaloc=etalocs[compo][trial]
            except:
                code_debugger()
            betaloc = beta[np.ix_(sidx, sidx)]
            if 'betalim' in kwargs:
                betalim=kwargs['betalim']
                betaloc=np.clip(betaloc,-betalim,betalim)

            betaerr=ifelse(len(bootdata['beta'])>1,0.001+np.std(bootdata['beta'].values,axis=0)[np.ix_(sidx, sidx)],np.ones(betaloc.shape))

            if  not 'bminfer' in sys.argv:
                # cidx=[species.index(c) for c in compo]
                try:
                    bminf, bvinf = hebbian_getbstats(allbeta[np.ix_(sidx,sidx)], etaloc, get_bm_from_beta=1,betaerr=betaerr)[:2]
                except:
                    code_debugger()
                bsinf,gaminf,tab=bvinf**.5,0,None
            else:
                bminf,bsinf,gaminf,tab=infer_bm(etaloc,S=Smax,table=tabs.get(compo,None),
                                         meanfield=kwargs.get('use_meanfield_bm', 'mf' in sys.argv))

            tabs[compo]=tab
            if kwargs.get('bm_per_replica', 'bminfer' in sys.argv ):
                bm=bminf
            else:
                bm=bmall
            if gamma is None:
                gam=gamall
            else:
                gam=gamma
            dic=correlcalc(etaloc,betaloc,gamma=gam,bm=bm,**kwargs)#,betaerr=betaerr,**kwargs)
            dic.update({'bmall':bmall,'bsall':bsall,'gamall':gamall,'bs':bsinf,'gamma':gaminf,'bmloc':bminf})
            if 'betalim' in kwargs:
                for k in dic:
                    if 'mean_emp' in k or 'mean_theo' in k:
                       dic[k]=np.clip(dic[k],-betalim,betalim)


            def indices(c):
                i = np.multiply.outer(range(len(c)), np.ones(len(c))).astype('int')
                return zip(offdiag(i), offdiag(i.T))

            sidx=np.array(sidx)
            dic.update({'composition':tuple(compo), 'richness':len(compo), 'eta':etaloc,'beta':betaloc})
            for z in ('corr_i','corr_j','corr_k','mean_i','mean_j'):
                dic[z]=sidx[list(dic[z])]
            results.append(dic)

    dfresults=pd.DataFrame(results)
    for result in results:
        dfm=pd.DataFrame({i:result[i] for i in result if 'mean' in i})
        dfc=pd.DataFrame({i:result[i] for i in result if 'corr' in i})

        adfm, adfc, runavgmat, theomat,rndmat = running_average(dfm, dfc, species,**kwargs)
        result.update({'runavgmat': runavgmat,'theomat':theomat,'mean_theo_runavg':adfm['mean_theo'].values, 'mean_emp_runavg':adfm['mean_emp'].values,
                       'corr_theo_runavg': adfc['corr_theo'].values, 'corr_emp_runavg': adfc['corr_emp'].values})

    dfresults=pd.DataFrame(results)
    if not len(dfresults):
        print 'BOOTDATA: NO RESULTS'
        code_debugger()
    return bootdata, dfresults



def show_bootstrap(bootdata,results,**kwargs):
    minusbm = kwargs.get('minusbm', 'bm' in sys.argv)
    scatter_alpha=0.5
    plotkw={}#'lowess':1,'ci':90}#'truncate':1, 'x_bins':8}
    import seaborn as sns
    species=bootdata['species'].values[0]#sorted(np.unique( np.concatenate(bootdata['eta'].values[0]['composition'].values ) ))
    S=len(species)
    path=Path(kwargs.get('path','.'))

    figdata = {}
    if not kwargs.get('redo',0):
        try:
            figdata=pickle.load(open(Path(path+kwargs.get('fname','results'))+'figdata.pickle','r'))
        except:
            pass
    figs={}
    theomat=None
    def perm(x):
        if 'noperm' in sys.argv:
            return x
        # Permutation of things to concatenate to break correlations induced by same bootstrap
        idx=range(len(x))
        return [x[i] for i in np.random.permutation(idx)]

    dfm=pd.DataFrame({x:np.concatenate(perm([np.array(z,dtype='float')-ifelse(minusbm,np.mean(nonan(np.array(z,dtype='float'))),0)
                    for z,b in results[[x,'bm']].values] )).astype('float') for x in ['mean_emp','mean_theo']})
    for x in results:
        if 'mean_' in x and not x in dfm:
            dfm[x] = np.concatenate(results[x].values)
    dfm['replica']= np.concatenate([np.ones(len(results['mean_emp'].values[i]),dtype='int' )*i for i in range(results.shape[0])])
    if 'bug' in sys.argv:
        # DEBUG FIGURE
        figs['debug']=plt.figure()
        plt.title(kwargs.get('title','')+' Check for bias due to emp/theo correlation across bootstraps'),plt.xlabel('Median(mean_emp)'),plt.ylabel('Median(mean_theo)')
        plt.scatter(*zip(*[ (np.median(nonan(np.array(x,dtype='float'))), np.median(nonan(np.array(y,dtype='float')))) for x,y in results[['mean_emp','mean_theo']].values  ]))

    # code_debugger()
    dfm=dfm[~np.isnan(np.sum(dfm.values,axis=1))]
    dfc=pd.DataFrame({x:np.concatenate(perm(results[x].values)).astype('float') for x in results if 'corr_' in x})
    dfc['replica']= np.concatenate([np.ones(len(results['corr_emp'].values[i]),dtype='int' )*i for i in range(results.shape[0])])
    dfm=dfm[~np.isnan(dfm['mean_emp'])].sort_values('mean_theo')
    dfc=dfc[~np.isnan(dfc['corr_emp'])].sort_values('corr_theo')
    
    #
    def sparsify(xs):
        bins=np.linspace(np.min(xs)-0.01,np.max(xs)+0.01,15)
        xbin=np.digitize(xs,bins)
        abund=np.array([np.sum(xbin==i) for i in range(len(bins)) ]) [xbin]
        idx=(np.random.random(xs.shape)<200./xs.shape[0])|(abund<10  )
        return idx



    fig, axes = plt.subplots(ncols=2,figsize=np.array(mpfig.rcParams['figure.figsize'])*(1.,.66)  )
    figs['reg']=fig
    mxs, mys = dfm['mean_theo'], dfm['mean_emp']
    mnidx  = sparsify(mxs)

    xs, ys = dfc['corr_theo'], dfc['corr_emp']
    idx = sparsify(xs)

    if kwargs.get('regplot',0):
        sns.regplot(ax=axes[0],x="mean_theo", y="mean_emp", data=dfm,scatter_kws={'alpha':scatter_alpha,'color':(.5,.7,.9)},line_kws={'color':'r','alpha':1} ,
                    **plotkw)
        #
        sns.regplot(ax=axes[1],x="corr_theo", y="corr_emp", data=dfc,scatter_kws={'alpha':scatter_alpha,'color':(.5,.7,.9)},line_kws={'color':'r','alpha':1},
                    **plotkw)

    else:
        axes[0].scatter(mxs[mnidx],mys[mnidx],**{'alpha':scatter_alpha,'color':(.5,.7,.9)})
        axes[1].scatter(xs[idx],ys[idx],**{'alpha':scatter_alpha,'color':(.5,.7,.9)})
        # print 'Starting savgol'
        # from scipy.signal import savgol_filter
        # axes[0].plot(mxs[mnidx],
        #     savgol_filter(mys[mnidx], 11, 3))

    print '======================== Binned running average ================================='

    beta = np.mean( bootdata['beta'].values, axis=0)
    betavar=np.var( bootdata['beta'].values, axis=0) # This is variance on the mean estimate of a given coefficient, not proper variance

    if 'adfm' in figdata:
        adfm,adfc,runavgmat,theomat,rndmat=[figdata[x] for x in ['adfm','adfc','runavgmat','theomat','rndmat']]
    else:
        thresh=max(.1, 10000./dfc.shape[0] )
        adfm,adfc,runavgmat,theomat,rndmat=running_average(dfm,dfc.iloc[np.where(np.random.random(dfc.shape[0])<thresh)[0] ],species,**kwargs)
        figdata.update({'adfm':adfm,'adfc':adfc,'runavgmat':runavgmat,'theomat':theomat,'rndmat':rndmat,'beta':beta, 'species':species})



    adfm=adfm.sort_values('mean_theo')
    adfc=adfc.sort_values('corr_theo')
    
    mnidx=sparsify(adfm['mean_theo'])
    axes[0].plot(adfm['mean_theo'].values[mnidx], adfm['mean_emp'].values[mnidx], alpha=.5, lw=4, color='k')
    axes[0].fill_between(adfm['mean_theo'].values[mnidx], (adfm['mean_emp']+adfm['mean_emp_05']).values[mnidx],
                         (adfm['mean_emp']+ adfm['mean_emp_95']).values[mnidx], alpha=.2,
                              lw=1, color='k')
    xs=adfc['corr_theo'].values
    cridx=sparsify(xs)
    axes[1].plot(xs[cridx], adfc['corr_emp'].values[cridx], alpha=.5, lw=4, color='k')
    axes[1].fill_between(xs[cridx],  (adfc['corr_emp']+adfc['corr_emp_05']).values[cridx],
                         (adfc['corr_emp']+adfc['corr_emp_95']).values[cridx], alpha=.2,
                              lw=1, color='k')
    if 'lowess' in sys.argv:
        axes[0].plot(*lowess(dfm['mean_emp'], dfm['mean_theo']).T, lw=5, color='g', alpha=.4)
        axes[1].plot(*lowess(dfc['corr_emp'], dfc['corr_theo']).T, lw=5, color='g', alpha=.4)


    ncompos=len(set([tuple(c) for c in results['composition'].values]))
    if ncompos > 1:
        raise Exception("show_bootstrap not functional when multiple compositions in results")
    try:
        betadiff = beta - runavgmat
    except:
        code_debugger()


    mergecomp=sorted(np.unique([x for c in results.composition.values for x in c])) 
    cidx=[species.index(c) for c in mergecomp]
    if  S>len(results['composition'].iloc[0]):
        betadiff=betadiff[np.ix_(cidx,cidx)]

    etameandf=pd.concat([pd.DataFrame([dict(list(zip(c,np.median(x,axis=1)))+[('composition',tuple(c)) ] )]) for x,c in results[['eta','composition']].values],ignore_index=True)
    etameandf=etameandf[mergecomp+['composition']]
    etamean=etameandf[list(mergecomp)].median().values

    xrg=dfm['mean_theo'].min(),dfm['mean_theo'].max()
    axes[0].plot(xrg,xrg,linestyle='--',color='r'),axes[0].set_title('Means')
    dxrg=xrg[1]-xrg[0]
    axes[0].set_ylim(ymin=xrg[0]-dxrg,ymax=xrg[1]+dxrg)
    xrg=dfc['corr_theo'].min(),dfc['corr_theo'].max()
    axes[1].plot(xrg,xrg,linestyle='--',color='r'),axes[1].set_title('Correlations')
    dxrg=xrg[1]-xrg[0]
    axes[1].set_ylim(ymin=xrg[0]-dxrg,ymax=xrg[1]+dxrg)

    if 'title' in kwargs:
        plt.suptitle(kwargs['title'])

    if kwargs.get('corrunavg', 0):
            """ TESTING COMPUTING CORRELATIONS AROUND RUNAVG INSTEAD OF BETATHEO =================================================  """

            Sb=betadiff.shape[0]
            diag = np.multiply.outer(np.ones(Sb), np.eye(Sb))

            def removeself(mat):
                S = mat.shape[0]
                ss = range(S)
                mat2 = [mat[i][np.ix_([s for s in ss if s != i], [s for s in ss if s != i])] for i in range(Sb)]
                return np.array(mat2)

            empirical = removeself(np.einsum('ij,ik->ijk', betadiff, betadiff))
            var = np.array([np.nanmean(empirical[i][np.eye(Sb - 1) != 0]) for i in range(Sb)]).reshape((-1, 1, 1))
            empirical /= var + 10 ** -15
            empirical -= removeself(diag)
            empirical=empirical.ravel()


            if ncompos==1:
                corrtheo = np.mean(results['corr_theo'], axis=0)
                prev = np.mean(results['corr_emp'], axis=0)
                newfig, newaxes = plt.subplots(ncols=2, figsize=np.array(mpfig.rcParams['figure.figsize']) * (1., .66))
                figs['corrunavg']=newfig
                sns.regplot(ax=newaxes[0], x=corrtheo,y=prev,
                            scatter_kws={'alpha': scatter_alpha, 'color': (.5, .7, .9)},
                            line_kws={'color': 'r', 'alpha': 1},
                            **plotkw),newaxes[0].set_title('Use beta - E[beta]_theo | slope: {:.2f}'.format(linregress(corrtheo,empirical)[0] ))
                sns.regplot(ax=newaxes[1], x=corrtheo,y=empirical,
                            scatter_kws={'alpha': scatter_alpha, 'color': (.5, .7, .9)},
                            line_kws={'color': 'r', 'alpha': 1},
                            **plotkw),newaxes[1].set_title('Use beta - runavg(beta) | slope: {:.2f}'.format(linregress(corrtheo,prev )[0]) )
                newaxes[0].plot(xrg, xrg, linestyle='--', color='k')
                newaxes[1].plot(xrg, xrg, linestyle='--', color='k')
            else:
                #OBSOLETE
                corrtheo = ( removeself(
                    - np.multiply.outer(1. / (np.sum(etamean ** 2) - etamean ** 2 + 0.0001), np.multiply.outer(etamean, etamean))) ).ravel()
                figs['corrunavg'] =plt.figure()
                sns.regplot(ax=plt.gca(), x=corrtheo, y=empirical,
                            scatter_kws={'alpha': scatter_alpha, 'color': (.5, .7, .9)},
                            line_kws={'color': 'r', 'alpha': 1},
                            **plotkw)
                plt.title('corr(beta_ij,beta_ik) using beta-runavg(beta)  | slope: {:.2f}'.format(linregress(corrtheo,empirical)[0] ))

            plt.figure(fig.number)


    figdata['mean']=results['mean'].values
    figdata['corr']=results['corr'].values
    figdata['composition']=results['composition'].values[0]

    if not theomat is None:
        """ COMPARE RUNAVG MATRIX WITH PREDICTED E[beta] """
        figs['matrix']=newfig=plt.figure()
        if 'title' in kwargs:
            plt.suptitle(kwargs['title'])


        from numpy.linalg import inv as lainv
        arg=np.argsort(etamean )
        def rgnot(x):
            r=range(beta.shape[0])
            r.remove(x)
            return r
        earg=[a for a in arg if etamean[a]>0 and (( beta[cidx[a],rgnot(cidx[a])] !=0).any() or ( beta[rgnot(cidx[a]),cidx[a]] !=0).any() )]
        arg=[cidx[a] for a in earg  ]
        aa=np.ix_(arg, arg)
        def sel(x):
            return nonan(offdiag(x[aa]))
        vmin,vmax=np.min(sel(np.minimum(runavgmat,theomat))), np.max(sel(np.maximum(runavgmat,theomat)))

        original_vmax=None

        bmm=results['bm'].median()
        if minusbm:
            vmax=max(np.abs([vmin,vmax]))
            vmin=-vmax
            cmap='seismic'
            original_vmax=vmax
            beta=beta-bmm
        else:
            #cmap='viridis'
            cmap='seismic'
            vmax=np.max(np.abs([bmm-vmin,vmax-bmm]))+bmm
            vmin=2*bmm-vmax
        np.fill_diagonal(beta,np.nan)
        np.fill_diagonal(runavgmat,np.nan)
        np.fill_diagonal(theomat,np.nan)
        np.fill_diagonal(rndmat,np.nan)
        if kwargs.get('remove_zeros',1):
            runavgmat[runavgmat==0]=np.nan
            theomat[theomat==0]=np.nan

        # arg=list(np.argsort([ np.mean(x[(~np.isnan(x))&(x>np.median(x)) ]) for x in theomat.T] ))
        #
        figdata[ 'beta_sort']=beta[aa]
        beta_std = np.std( bootdata['beta'].values, axis=0)
        figdata[ 'beta_sort_std']=beta_std[aa]
        figdata[ 'eta_sort']=etamean[earg]
        figdata[ 'runavgmat_sort']=runavgmat[aa]
        figdata[ 'rndmat_sort']=rndmat[aa]
        figdata[ 'theomat_sort']=theomat[aa]
        figdata['species_sort']=[species[s] for s in arg]

        # code_debugger()

        plt.subplot(331)
        betamax=np.max(nonan(offdiag(np.abs(beta[aa]))))
        vmin,vmax=2*(bmm) -betamax,betamax
        plt.imshow(beta[aa],cmap=cmap,vmin=2*bmm -betamax,vmax=betamax),plt.colorbar(),plt.title("Measured"),plt.yticks(range(len(arg)),[species[s][:3] for s in arg]),plt.xticks([])
        plt.subplot(332)
        plt.imshow(runavgmat[aa],vmin=vmin,vmax=vmax,cmap=cmap),plt.colorbar(),plt.title("Run avg"),plt.xticks([]),plt.yticks([])
        plt.subplot(333)
        plt.imshow(theomat[aa],vmin=vmin,vmax=vmax,cmap=cmap),plt.colorbar(),plt.title("Predicted"),plt.xticks([]),plt.yticks([])
        plt.subplot(337)
        b=beta[aa]
        sec=2
        bsplit=[[np.median(nonan(x)) for x in np.array_split(a,sec,axis=1) ] for a in np.array_split(b,sec,axis=0)]
        plt.imshow(bsplit,vmin=vmin,vmax=vmax,cmap=cmap),plt.colorbar(),plt.title("Quadrant avg"),plt.xticks([]),plt.yticks([])
        # plt.imshow(rndmat[aa],vmin=vmin,vmax=vmax,cmap=cmap),plt.colorbar(),plt.title("Randomized"),plt.xticks([]),plt.yticks([])

        # print 'empirical,runavg,theo',np.mean(offdiag(beta[aa])), np.mean(offdiag(runavgmat[aa])),np.mean(offdiag(theomat[aa]))
        plt.subplot(338)
        from scipy.ndimage import convolve
        b=beta[np.ix_(arg,arg)].copy()
        bmm=results['bmall'].median()
        np.fill_diagonal(b,bmm)
        b=convolve(b,np.ones((5,5))/25.,mode='constant',cval=bmm)
        np.fill_diagonal(b,np.nan)
        plt.imshow(b,vmin=vmin,vmax=vmax,cmap=cmap),plt.colorbar(),plt.title("Naive smooth"),plt.xticks([]),plt.yticks([])

        diffs=[(beta-theomat,'Meas - Pred'),(beta-runavgmat,'Meas - Runavg'),
               (runavgmat-theomat,'Runavg - Pred'), (np.median(list(results['runavgmat'].values),axis=0),'Median boot'),
               ]
        for z,diff in enumerate(diffs[:-1] ):
            diff,titl=diff
            diff=diff[aa]
            plt.subplot(334+z)
            if minusbm and 'Runavg - Pred' in titl:
                vmax=original_vmax
            else:
                vmax=np.max(np.abs(nonan(diff)))
            plt.imshow(diff,cmap='seismic',vmin=-vmax,vmax=vmax),plt.colorbar(),plt.title(titl),plt.xticks([]),plt.yticks([])

        plt.subplot(339)
        theofill=theomat[aa]
        dia=list((np.diag(theofill,1)+np.diag(theofill,-1))/2)
        dia.append(dia[-1])
        theofill[np.eye(theofill.shape[0])>0]=dia
        vmin,vmax=2*(bmm) -betamax,betamax

        plt.colorbar(plt.imshow(beta[aa],cmap=cmap,vmin=vmin,vmax=vmax)),plt.title("Iso curves"),plt.xticks(range(len(arg)),[species[s][:3] for s in arg],rotation=90),plt.yticks([])
        # levels=np.concatenate([np.linspace(np.min(theofill),0,theofill.shape[0]),np.linspace(0,np.max(theofill),theofill.shape[0])[1:] ])
        levels=np.unique(sorted(theofill.ravel()))
        levels=levels[::len(levels)/8]

        from scipy.ndimage.filters import gaussian_filter
        plt.contour(gaussian_filter(theofill,1), linewidths=8, colors='w',alpha=.5, linestyles='solid', levels=levels )
        plt.contour(gaussian_filter(theofill,1), linewidths=4, colors='k',alpha=.5, linestyles='solid', levels=levels )

        cidx=[species.index(c) for c in mergecomp]
        Kmono = np.mean([z for z in bootdata['Kmono'].values], axis=0)[ cidx]
        Nmean=etameandf.copy()
        Nmean[list(mergecomp)]=etameandf[list(mergecomp)]*Kmono

        ccfigs,ccmeas=cc_analysis(bootdata, mergecomp, Nmean, Kmono, etameandf,beta=beta[np.ix_(cidx,cidx)],runavgmat=runavgmat[np.ix_(cidx,cidx)],
                           betatheo=theomat[np.ix_(cidx,cidx)],bm=results['bmall'].median(),bs=results['bsall'].median() ,betavar=betavar[np.ix_(cidx,cidx)],all_plots=0,
                                  title=kwargs['title'])
        for i,j in ccfigs.items():
            figs[i]=j
        figdata['slopes_per_species']=ccmeas


    rpath=Path(path+kwargs.get('fname','results'))
    rpath.mkdir()
    if kwargs.get('save',None):
        for fig in figs:
            figs[fig].savefig(rpath+fig+'.svg' )

    if not kwargs.get('hold',0):
        plt.show()

    if kwargs.get('debug',0):
        code_debugger()

    pickle.dump(figdata,open(rpath+'figdata.pickle','w'))
    return figdata

def bilin(etai,etaj,mat):
    def bi(x,A,B,C,D):
        return A*x[0]*x[1] + B*x[0] + C*x[1] +D
    res=curve_fit(bi,np.array([etai,etaj]),mat )
    return res[0][0]



def exp_simus(exp,results,species,path='.',runbins=5):
    """Simulations with same etas and statistics as experiment"""
    sim={}
    basepath=Path(path)
    etamean = np.mean([np.mean(e,axis=1) for e in results['eta'].values], axis=0)
    S=len(etamean)
    bm,bs=results[['bmall','bsall']].mean()
    for mode in ['gauss','cc']: #'worst',
        path=basepath+Path('compar_{}'.format(mode))
        path.mkdir()
        if 'bootdata.json' in os.listdir(path) and not 'resim' in sys.argv and not 'rerun' in sys.argv :
            bootdata, results = pd.read_json(path + 'bootdata.json'), pd.read_json(path + 'bootresults.json')
            for d in (bootdata, results):
                for key in d:
                    if key in ('composition', 'species'):
                        d[key] = [tuple(c) for c in d[key]]
                        continue
                    test = d[key].values[0]
                    if hasattr(test, 'keys'):
                        pass
                    elif hasattr(test, '__iter__') and not (
                                isinstance(test, tuple) or isinstance(test, basestring)):
                        d[key] = [np.array(m, dtype='float') for m in d[key].values]
        else:
            bootdata, results =simus(betamode=mode, etamode='use',use_eta=etamean, bm=bm,bs=bs,
                                  ncompo=1, niter=10, S=S, runbins=runbins,death=0)
            bootdata.to_json(path + 'bootdata.json'), results.to_json(path + 'bootresults.json')
        results['composition'] = [tuple(species) for c in results.composition]
        if len(results['composition'].values[0]) != len(results['eta'].values[0]):
            code_debugger()
        sim[mode]=(bootdata,results)
    return sim


def pattern_score(experiments,trials=10,path='.',**kwargs):
    import seaborn as sns
    from scipy.linalg import norm as lanorm
    from scipy.special import erf
    from scipy.stats import spearmanr

    alldats={'base':[], 'replica':[],'worst':[] ,'gauss':[],
             'cc':[]}
    figs={}



    for exp in experiments:
        df,bootdata,results,figdatas,expsimus=experiments[exp]
        testcompos = sorted(set([tuple(c) for c in results['composition'].values]))
        #
        species=bootdata['species'].values[0]
        # FOR RUNAVG PER BOOTSTRAPPED REPLICA

        betaerr = np.nanstd(list(bootdata['beta'].values), axis=0)
        betaerr[np.isnan(betaerr)] = np.max(betaerr[~np.isnan(betaerr)])
        betaerr[betaerr == 0] = np.max(betaerr)

        # FOR RUNAVG OVER ALL REPLICAS
        for i,figdata in enumerate(figdatas):
            ifig=i
            bs = results['bsall'].mean()
            compo=figdata.get('composition',testcompos[0])
            cidx=[species.index(ss) for ss in compo]
            locbetaerr=betaerr[np.ix_(cidx,cidx)]
            comparisons = {'replica': results[results.composition==compo]}
            if hasattr(expsimus,'keys'):
                for s in expsimus:
                    comparisons[s] = expsimus[s][1]
            else:
                for s in expsimus[i]:
                    comparisons[s] = expsimus[i][s][1]
            for s in comparisons:
                comparisons[s]=comparisons[s][comparisons[s].composition==compo]

            """ MEANS ================ """

            runavg, theomat, beta, rndmat, etamean, spsort = (figdata['runavgmat_sort'], figdata['theomat_sort'],
                                                             figdata['beta_sort'], figdata['rndmat_sort'], figdata[
                                                                 'eta_sort'], figdata['species_sort'])
            sidx = [compo.index(s) for s in spsort if s in compo]
            sort = np.ix_(sidx,sidx)
            uidx = [spsort.index(s) for s in compo]
            unsort = np.ix_(uidx, uidx)

            corrempsuffix=['','_runavg'][0]
            print '      NB: for correlations, comparing theo to emp{}'.format(corrempsuffix)
            theomat=theomat.copy()
            np.fill_diagonal(theomat,0)

            cemp=np.nanmean(list(comparisons['replica']['corr_emp'+corrempsuffix].values), axis=0)
            ctheo, corri, corrj, corrk=[np.nanmean(list(comparisons['replica'][z].values), axis=0) for z in ('corr_theo','corr_i','corr_j','corr_k')]


            # adfc = figdata['adfc'].groupby(['corr_i', 'corr_j', 'corr_k']).median().reset_index()

            tmp=pd.DataFrame({ii:[jj] for ii,jj in {'beta':beta[unsort], 'runavgmat':runavg[unsort], 'theomat':theomat[unsort],
                              'corr_emp'+corrempsuffix:cemp, 'corr_theo':ctheo,'corr_i':corri,'corr_j':corrj, 'corr_k':corrk,
                 'eta':etamean[uidx].reshape(etamean.shape+(1,)),'bsall':comparisons['replica']['bsall'].mean() }.items()} )
            comparisons['base']=tmp
            # z=tmp['beta'].values[0]
            # plt.close('all')
            # [plt.scatter(offdiag(z), offdiag(w)) for w in comparisons['replica']['beta'].values]
            # plt.plot(offdiag(z),np.median([offdiag(w) for w in comparisons['replica']['beta'].values],axis=0)  )
            # plt.show()

            for loclabel, locresults in comparisons.items():
                bs = locresults['bsall'].mean()
                if  'output' in sys.argv:
                    opath=Path(path)+Path('output')
                    opath.mkdir()
                    # print opath
                    compnb = ''
                    if len(figdatas) > 1:
                        compnb = '_{}'.format(ifig)
                    if 'base' in loclabel:
                        beta=np.median(locresults['beta'].values,axis=0)
                        np.savetxt(opath+'{}{}_matrix.csv'.format(exp,compnb),beta )
                        np.savetxt(opath+'{}{}_matrix_err.csv'.format(exp,compnb),betaerr )
                        np.savetxt(opath+'{}{}_equilibrium.csv'.format(exp,compnb),etamean )
                    if 'gauss' in loclabel or 'cc' in loclabel:
                        fname=opath+'{}{}_{}.csv'.format(exp,compnb,loclabel)
                        fil=open(fname,'w')
                        fil.close()
                        for imat,mat in enumerate(locresults['beta'].values ):
                            if 'cc' in loclabel:
                                mat = snap_beta(mat, np.median(locresults['eta'].values[imat],axis=1) * np.max(etamean) )
                            # print loclabel, np.dot(mat,etamean)
                            fil=open(fname,'a')
                            np.savetxt(fil,mat )
                            fil.write('\n######################\n\n')
                            fil.close()


                for beta, runavg, theomat, cemp, ctheo, corri,corrj,corrk, etamean in locresults[
                    ['beta', 'runavgmat', 'theomat', 'corr_emp'+corrempsuffix, 'corr_theo','corr_i','corr_j','corr_k',
                     'eta']].values:

                    """NB: SORT BY INCREASING ETA\\
                    All the variances come from cc_analysis where species are sorted by abundance"""



                    cemp = np.array(cemp)
                    etamean = np.median(etamean, axis=1)

                    if 'snap' in sys.argv or 'cc' in loclabel:
                        beta=snap_beta(beta,etamean)
                        np.fill_diagonal(beta,np.nan)
                    elif 'cc' in loclabel:
                        etamean=etamean/np.max(etamean)

                    a = [species.index(s) for s in compo]
                    try:
                        assert len(a) == len(etamean)
                        # raise
                    except:
                        code_debugger()
                    a = np.ix_(a, a)

                    beta,runavg,theomat=[mm[sort] for mm in (beta,runavg,theomat)]
                    etamean = etamean[sidx]

                    etai = np.multiply.outer(etamean, np.ones(len(etamean)))

                    offtheo = offdiag(theomat)
                    offrun = offdiag(runavg)
                    if kwargs.get('remove_zeros', not 'simu' in exp and not 'cc' in loclabel):
                        rem0 = (offtheo != 0) & (~np.isnan(offtheo))
                    else:
                        rem0 = (~np.isnan(offtheo))
                    rem0=rem0& (~np.isnan(offrun))& (~np.isnan(offdiag(beta)))
                    if kwargs.get('remove_err', 1):
                        rem0 = rem0 & (offdiag(locbetaerr[sort]) < bs )

                    mask = (np.eye(beta.shape[0]) == 0)
                    mask[mask==True] = rem0

                    if np.isnan(offdiag(theomat)).any():
                        print "NAN IN THEOMAT", exp, theomat
                        #code_debugger()
                    x, y = offrun[rem0], offtheo[rem0]
                    if len(figdatas) > 1 and not 'simu' in exp:
                        suffix = '_{}'.format(testcompos.index(compo))
                    else:
                        suffix = ''
                    sij = np.arange(len(compo))
                    i = np.multiply.outer(sij, np.ones(len(compo)))
                    i, j = offdiag(i)[rem0], offdiag(i.T)[rem0]
                    #
                    x, xalt, y = offdiag(runavg)[rem0], offdiag(beta)[rem0], offdiag(theomat)[rem0]
                    etai, etaj = offdiag(etai)[rem0], offdiag(etai.T)[rem0]
                    try:
                        slopx, slopy = bilin(etai, etaj, xalt), bilin(etai, etaj, y)
                    except:
                        code_debugger()
                    corri,corrj,corrk=[z.ravel() for z in (corri,corrj,corrk)]

                    if loclabel in ['base','replica']:
                        corri,corrj,corrk=[  np.array( [ compo.index(species[int(np.round(z))]) for z in corrc] ) for corrc in (corri,corrj,corrk) ]

                    try:
                        cetai=etamean[list(int(z) for z in corri )]
                        cetaj=etamean[list(int(z) for z in corrj )]
                        cetak=etamean[list(int(z) for z in corrk )]
                    except:
                        code_debugger()
                    corrjk=corrj*corrk
                    cetajk=cetaj*cetak

                    doboot=('base' in loclabel) and ( not 'simu' in exp) 
                    cgood=(~np.isnan(cemp))
                    ctheo,cemp,cetai,cetajk,corri,corrjk=[z[cgood] for z in (ctheo,cemp,cetai,cetajk,corri,corrjk)]
                    for trial in range(ifelse(doboot, 10, 1)):
                        dic = {'mean_pearson': np.corrcoef(x, y)[0, 1],
                               'mean_spearman': spearmanr(x, y)[0],
                               'mean_reldist': reldist(x, y, typ='', boot=doboot),
                               'corr_pearson': np.corrcoef(ctheo, cemp)[0, 1],
                               'corr_spearman': spearmanr(ctheo, cemp)[0],
                               'corr_reldist': reldist(ctheo, cemp, typ='', boot=doboot),
                               'corr_slopei': slopedist(ctheo,cemp,cetai,cetai,per_species=(corrjk,corrjk, sij)),
                               'corr_slopejk': slopedist(ctheo,cemp,cetajk,cetajk,per_species=(corri,corri, sij)),
                               'mean_sloperow': slopedist(x, y, etai, etaj, boot=doboot, per_species=(i, j, sij)),
                               'mean_slopecol': slopedist(x, y, etaj, etai, boot=doboot, per_species=(j, i, sij)),
                               'mean_sloperow_beta': slopedist(xalt, y, etai, etaj, boot=doboot, per_species=(i, j, sij)),
                               'mean_slopecol_beta': slopedist(xalt, y, etaj, etai, boot=doboot, per_species=(j, i, sij)),
                               'mean_meanrow_beta': rowdist(xalt, y, i,species=sij),
                               'mean_meancol_beta': rowdist(xalt, y, j,species=sij),
                               'mean_bilin':slopx,# 1 - 2 * np.abs(slopx - slopy) / (np.abs(slopx) + np.abs(slopy)),
                               'mean_smooth': smoothness(beta,mask=mask, bs=bs),
                               'mean_mn': np.mean(xalt),
                               'mean_sd': np.std(xalt),

                               }
                        # dic['mean_mindist']=np.min([dic['mean_'+z] for z in ['reldist','sloperow','slopecol']],axis=0)
                        if np.isnan(dic.values()).any():
                            code_debugger()
                        dic['exp']= exp + suffix

                    try:

                        dic['mean_loglike_row'] = slopedist(xalt, y, etai, etaj, typ='loglike',
                                var=figdata['slopes_per_species']['rowslope']['var'] ,
                                boot=0, per_species=(i, j, sij))  # /figdata['slopes_per_species']['colslope']['reldist_gen50']
                        dic['mean_loglike_col'] = slopedist(xalt, y, etaj, etai, typ='loglike',
                                var=figdata['slopes_per_species']['colslope']['var'],
                                boot=0, per_species=(j, i, sij))  # /figdata['slopes_per_species']['rowslope']['reldist_gen50']


                        dic['mean_loglike_row_mn'] = rowdist(xalt, y, i, typ='loglike',
                                var=figdata['slopes_per_species']['rowmean']['var'] ,
                                boot=0,species=sij)  # /figdata['slopes_per_species']['colslope']['reldist_gen50']
                        dic['mean_loglike_col_mn'] = rowdist(xalt, y, j, typ='loglike',
                                var=figdata['slopes_per_species']['colmean']['var'],
                                boot=0,species=sij)  # /figdata['slopes_per_species']['rowslope']['reldist_gen50']
                    except:
                        print "ERROR",loclabel
                        code_debugger()

                    varb=figdata['slopes_per_species']['overall']['var']
                    if len(varb.shape)==1:
                        varb=varb[rem0]
                    else:
                        varb=varb[np.ix_(rem0.ravel(), rem0.ravel())]
                    try:
                        loglike=reldist(xalt,y,boot= 0,strip=0, typ='loglike',var= varb )
                        dic['mean_loglike'] = loglike-figdata['slopes_per_species']['overall']['loglikesimu']
                    except:
                        code_debugger()

                    # if 'base' in loclabel:
                    #     code_debugger()

                    alldats[loclabel].append(dic)

            # code_debugger()

    for k in alldats:
        if not len(alldats[k]):
            continue
        dat=pd.DataFrame(alldats[k])
        dat['exp'] = [z.split('_')[0].split('+')[0] if kwargs.get('merge', 0) else z for z in
                      dat['exp'].values]  # Merge various compositions and simus
        alldats[k] =dat

    #
    dat=alldats['base']
    if not dat.shape[0]:
        return
    replica_dat=alldats['replica']
    cc_dat=alldats['cc']
    gauss_dat=alldats['gauss']
    # dat['expscore']=np.concatenate([[ '{} (p={:.3g})'.format(exp,gp['pval'].max() ) for z in range(gp.shape[0])] for exp,gp in dat.groupby('exp')] )
    #


    whis=[2.5,97.5]



    from matplotlib.gridspec import GridSpec
    with sns.axes_style("white"):
        sns.set_palette('deep')

        toshows= [('full', {'mean': ('bilin',
                   'sloperow_beta','slopecol_beta',
                  'reldist',
                  'loglike_col', 'loglike_row',
                  'loglike_col_mn', 'loglike_row_mn',
                  'loglike', 'smooth'),
         'corr': ('reldist','slopei','slopejk')}), ]
         
        if kwargs.get('articleplots'):
             toshows+=[ ('show', {'corr':['reldist'], 'mean':['bilin','smooth',],'merge':True }) ]


        showdic = {'mean': '', 'sloperow_beta': 'row slope', 'slopecol_beta': 'col slope', 'bilin': 'bilinear',
               'loglike_col': 'column slope LL', 'loglike_row': 'row slope LL',
               'loglike_col_mn': 'column mean LL', 'loglike_row_mn': 'row mean LL', 'loglike': 'overall LL',
               'smooth': 'smoothness',
               'slopei': 'slope with $\eta_i$', 'slopejk': 'slope with $\eta_{j}\eta_k$'}
        ###  ==================== TABLE OF SCORES
        total=None
        for name,toshow in toshows:
            # print toshow
            locdat = alldats['base'].copy().sort_values('exp')
            locreplica_dat = pd.concat(
            [replica_dat[[exp in e for e in replica_dat['exp']]].copy() for exp in locdat['exp'].unique()])

            for supti, heredat in [('Final matrix', locdat),]:
                                  # ('Per bootstrapped replica', locreplica_dat)]:
                # gpbe = sorted(gpbe, key=lambda e: e[1]['exppos'].values[0])

                heredat['exp'] = [z.split('_')[0].split('+')[0] if toshow.get('merge', 0) else z for z in
                              heredat['exp'].values]  # Merge various compositions and simus

                # hlist=['mean_{}'.format(zz) for zz in toshows[0][1]['mean']]+['corr_{}'.format(zz) for zz in toshows[0][1]['corr']]
                # heredat['mean_total']=[np.mean(hrow>0) for hrow in heredat[hlist].values]
                gpbe = heredat.groupby('exp')
                loctab = []
                for exp, gp in gpbe:
                    # suffixes=[]
                    for prefix in ['mean','corr']:
                        for suffix in toshow[prefix]:
                            xmlab = '{}_'.format(prefix) + suffix
                            locval = gp[xmlab].mean()
                            # if xmlab =='mean_total':
                            #     score=locval
                            #     loctab.append({'exp': exp, 'metric': 'total', 'basescore': score,
                            #                    'good': score, 'score': score,
                            #                    'discrimination': 1,})
                            #     code_debugger()
                            #     continue
                            locsd = gp[xmlab].std()
                            if np.isnan(locsd):
                                locsd = 0
                            worst = alldats['cc'].copy()
                            compar = alldats['gauss'].copy()
                            worst['exp'] = [z.split('_')[0].split('+')[0] if toshow.get('merge', 0) else z for z in
                                            worst['exp'].values]
                            compar['exp'] = [z.split('_')[0].split('+')[0] if toshow.get('merge', 0) else z for z in
                                             compar['exp'].values]

                            worst = worst[worst.exp == exp]
                            compar = compar[compar.exp == exp]
                            refval, worstval = compar[xmlab].mean(), worst[xmlab].mean()
                            refsd, worstsd = compar[xmlab].std(), worst[xmlab].std()

                            expo=1.
                            score = np.abs(locval - refval) / (2. * (refsd**expo + 0*locsd**expo)**(1./expo)  )  # COMPARISON to 2 sigma #  np.abs(refval-worstval)
                            discrimination = 1 / (1 + max(  (worstsd**expo + refsd**expo )**(1./expo) / np.abs(worstval - refval),
                                                            (refsd ** expo + locsd ** expo) ** (
                                                            1. / expo) / np.abs(worstval - refval)) )

                            logistic_score= 1- 2*(score>1.) #  1- 2/(1+ 1./score )
                            # logistic_discrimination= 1/(1+ 2*(worstsd+refsd)/np.abs(worstval - refval)   )
                            # if 'jk' in suffix:
                            #     code_debugger()

                            # scoredisc=np.log(score)*discrimination
                            loctab.append({'exp': exp, 'metric': prefix+' '+suffix,'basescore':score,'good':logistic_score, 'score': logistic_score*discrimination,
                                           'discrimination': discrimination})
                loctab=pd.DataFrame(loctab)
                suffixes=sorted(set(loctab.metric))
                # code_debugger()
                scores = loctab.pivot_table(index='exp', columns='metric')['score']
                # scores = scores[[s for s in suffixes if not True in [z in s  for z in noscore]]]
                if name == 'full':
                    total = np.mean(scores>0,axis=1)
                    scores = scores.iloc[np.argsort([np.mean(np.clip(scores.loc[ee].values, -2, 2)) for ee in scores.index])]
                elif name == 'show':
                    scores['total']=total
                    scores = scores.iloc[np.argsort(scores.total)]

                # scores=scores[[ 'bilin','sloperow_beta','slopecol_beta','meanrow_beta', 'meancol_beta','loglike','smooth' ]]
                figs['score_table_'+name+' '+supti] = plt.figure(
                    figsize=np.array(mpfig.rcParams['figure.figsize']) * (len(suffixes) / 2., 1.))
                plt.suptitle(supti+ ' '+name)
                # ax1 = plt.subplot(211)
                # ax2 = plt.subplot(212)
                ax1 = plt.gca()
                vmax = min(2., np.max(np.abs(scores.values)))
                sns.heatmap(scores.values, annot=True, fmt=".2g", linewidths=.5, ax=ax1, vmin=-vmax, vmax=vmax,
                            cmap='seismic_r')
                # code_debugger()
                plt.yticks(ax1.get_yticks(), scores.index, rotation=0)
                plt.xticks(ax1.get_xticks(),[' '.join([showdic.get(kk,kk) for kk in k.split(' ') ]) for k in scores.keys()], rotation=45)
                # ax1.set_title("Score")
        
        ##### ========================== BOXPLOTS
        toshow={'mean': ('bilin',
                   'sloperow_beta','slopecol_beta',
                     # 'sloperow','slopecol',
                  # 'meanrow_beta', 'meancol_beta',
                  'reldist',
                  'loglike_col', 'loglike_row',
                  'loglike_col_mn', 'loglike_row_mn',
                  'loglike', 'smooth',
                         # 'mn', 'sd',
                         ),
         'corr': ('reldist','slopei','slopejk')}

        for prefix,suffixes in toshow.items():
            if kwargs.get("replicaplots"):
                figs['{}_scores_replica'.format(prefix) ]=figrep=plt.figure(figsize=np.array(mpfig.rcParams['figure.figsize']) * (len(suffixes)/2., 1.))
                plt.suptitle('{} (per bootstrapped replica)'.format(prefix) )
            figs['{}_scores'.format(prefix) ]=fig=plt.figure(figsize=np.array(mpfig.rcParams['figure.figsize']) * (len(suffixes)/2., 1.))
            plt.suptitle('{} (final matrix)'.format(prefix) )
            spec = GridSpec(1, len(suffixes))

            xmlab = '{}_'.format(prefix) + suffixes[0]
            locdat = dat.sort_values('exp')
            locdat['exppos'] = np.concatenate(
                [[ifelse('simu' in exp, gp[xmlab].median(), gp[xmlab].median()) for z in range(gp.shape[0])] for exp, gp
                 in locdat.groupby('exp')])
            # dat['exppos']=np.concatenate([[gp[xmlab].median() for z in range(gp.shape[0])] for exp,gp in dat.groupby('exp')] )
            locdat = locdat.sort_values('exppos')
            locreplica_dat = pd.concat(
                [replica_dat[[exp in e for e in replica_dat['exp']]] for exp in locdat['exp'].unique()])



            ### FULL COMPARISON
            for idx,suffix in enumerate(suffixes):
                locfigs=[(locdat, fig)]
                xmlab='{}_'.format(prefix)+suffix
                if xmlab in locreplica_dat and kwargs.get("replicaplots"):
                    locfigs.append((locreplica_dat,figrep))

                for dathere,f in locfigs:
                    plt.figure(f.number)
                    ax1=plt.subplot(spec[idx])
                    # from brokenaxes import brokenaxes
                    # bax = brokenaxes(xlims=((-1, -7.), (.7, 1.)), subplot_spec=ax1)
                    sns.boxplot(data=dathere, y='exp', x=xmlab, whis=whis, ax=ax1)
                    plt.xlabel(' '.join([showdic.get(kk,kk) for kk in suffix.split(' ') ]) ,rotation=45)
                    if xmlab in cc_dat.keys():
                        for pos, text in zip(ax1.get_yticks(),ax1.get_yticklabels( )):
                            for datsrc, color in [('cc','r'),('gauss','b')]: #('worst','r'),
                                locexp=text.get_text()
                                locworst=alldats[datsrc][alldats[datsrc].exp==locexp]
                                ax1.errorbar(locworst[xmlab].mean(),pos,xerr=locworst[xmlab].std(),color=color, fmt='o',capsize=15,)
                            ax1.scatter(dathere.loc[dathere.exp==locexp,xmlab].mean(), pos,s=100 , c='k', )

                    if idx > 0:
                        plt.yticks([])
                    if  True in [ss in suffix for ss in ['reldist',] ]:#[ 'loglike', 'smooth', 'mean', 'mn' ,'sd','']]:
                        ax1.set_xlim(xmin=-1.1, xmax=1.1)



    tmp=(np.random.random(locdat['mean_sloperow'].values.shape)<.1)
    if len(tmp)>5:
        figs['metrics']=plt.figure()
        plt.suptitle('Agreement between metrics')
        plt.subplot(121),plt.title('sloperow  vs  reldist')
        scatter(locdat['mean_sloperow'].values[tmp],locdat['mean_reldist'].values[tmp],alpha=.1,hold=1,c='g'),plt.axhline(0),plt.axvline(0),plot([-1,1],[-1,1],hold=1)
        plt.subplot(122),plt.title('sloperow  vs  slopecol')
        scatter(locdat['mean_sloperow'].values[tmp],locdat['mean_slopecol'].values[tmp],alpha=.1,hold=1,c='g'),plt.axhline(0),plt.axvline(0),plot([-1,1],[-1,1],hold=1)

    if kwargs.get('debug',0):
        code_debugger()

    rpath=Path(path+kwargs.get('fname','results'))
    rpath.mkdir()
    if kwargs.get('save',None):
        for fig in figs:
            figs[fig].savefig(rpath+fig+'.svg' )

    return dat




def simus(ncompo=1,niter=1,S=16,bm=.2,bs=.15,betamode='',etamode='cc',gamma=0,noise=0,death=0.0,**kwargs):
    import numpy.linalg as la

    print "SIMU",etamode,betamode
    results=[]
    allS=S*ncompo
    betas = np.zeros((allS,allS))
    allspecies=[str(x) for x in range(allS)]
    if allS>1000:
        raise Exception('simus: matrix is too large')
    for comp in range(ncompo):
        ccmat = make_cc(S)
        if etamode=='use':
            eta=kwargs.get('use_eta',None)
            if eta is None:
                raise Exception("Missing etas in simus")
        elif etamode=='cc':
            eta = np.dot(la.inv(ccmat), np.ones(S))
        elif etamode=='pow':
            eta=0.9**np.arange(S)
        elif etamode=='simu':
            b=np.random.normal(bm,bs,(S,S))
            np.fill_diagonal(b,0)
            eta = dosimu(-b, np.ones(S), np.ones(S), noise=0,tmax=100)[-1][-1]
        elif etamode=='uniform':
            eta=np.random.uniform(0,1,S)
        else:
            eta=np.random.normal(1.,.3,S)/4
        alive=(eta>death)
        eta=eta[alive]
        ccmat=ccmat[np.ix_(alive,alive)]
        S=len(eta)
        ccmat=ccmat[np.ix_(np.argsort(eta),np.argsort(eta))]
        eta=np.sort(eta)
        bm_MF = (1. / np.mean(eta) - 1) / S * 1.1
        if etamode!='use':
            bm = hebbian_stablebm(eta)
            if etamode!='simu':
                bs=np.sqrt( (1- np.mean(eta)**2/np.mean(eta**2))/S )
                print 'BM {} BS {}'.format(bm,bs)


        species = [str(x+S*comp) for x in range(S)]

        trial=0

        while trial < niter:
            etaloc=eta + ifelse(noise,np.random.normal(0,noise,eta.shape),0)
            alive=np.where(etaloc>death)[0]
            compo=[species[i] for i in alive]
            etaloc=etaloc[etaloc>0]
            if betamode=='ccpert':
                betaloc=ccmat[np.ix_(alive,alive)]+create_synthetic_network_LV_gam0(etaloc, 0, bs,val=0)
                np.fill_diagonal(betaloc,0)
                betaloc=(betaloc.T*(1-etaloc)/(1-2*etaloc)).T+np.eye(len(alive))
            elif betamode=='gauss':
                betaloc = create_synthetic_network_LV_gam0(etaloc, bm, bs)
            elif betamode=='gausshigh':
                betaloc = create_synthetic_network_LV_gam0(etaloc, bm*2., bs)
            elif betamode=='worst':
                Slive=len(alive)
                epsilon=0.05
                betatheo = hebbian_getbstats(np.ones((Slive,Slive)), etaloc)[-1]
                betaloc=create_synthetic_network_LV_gam0(etaloc, bm, bs)
                worst=make_worst(etaloc,bm,bs)
                betaloc=worst

            elif betamode=='cc':
                etaloc=etaloc/np.max(etaloc)
                print trial, comp
                try:
                    betaloc=make_cctradeoff(etaloc)
                    if betaloc is None:
                        print 'FAILED'
                        betaloc=np.random.random((len(etaloc),len(etaloc)))
                        np.fill_diagonal(betaloc,1)
                    else:
                        betaloc=betaloc[-1]
                except:
                    continue
            trial+=1
            if betaloc is None:
                continue
            sidx=[allspecies.index(c) for c in compo]
            betas[np.ix_(sidx,sidx)]+=betaloc/niter
            dic=correlcalc(etaloc,betaloc,gamma=gamma,remove_zeros=0,**kwargs)

            def indices(c):
                i = np.multiply.outer(range(len(c)), np.ones(len(c))).astype('int')
                return zip(offdiag(i), offdiag(i.T))

            sidx=np.array(sidx)
            for z in ('corr_i','corr_j','corr_k','mean_i','mean_j'):
                dic[z]=sidx[list(dic[z])]

            dic.update({'composition':tuple(compo), 'richness':len(compo), 'eta':etaloc.reshape(-1,1),'beta':betaloc })
            dic['bmall']=bm
            dic['bsall']=bs
            results.append(dic)


    for result in results:
        dfm=pd.DataFrame({i:result[i] for i in result if 'mean' in i})
        dfc=pd.DataFrame({i:result[i] for i in result if 'corr' in i})
        adfm, adfc, runavgmat, theomat,rndmat = running_average(dfm, dfc, allspecies,**kwargs)
        result.update({'runavgmat': runavgmat,'theomat':theomat,'mean_theo_runavg':adfm['mean_theo'].values, 'mean_emp_runavg':adfm['mean_emp'].values,
                       'corr_theo_runavg': adfc['corr_theo'].values, 'corr_emp_runavg': adfc['corr_emp'].values})
        # code_debugger()
    dfresults=pd.DataFrame(results)

    bootdata=pd.DataFrame({'beta':[betas],'species':[allspecies],'Kmono':[np.ones(len(allspecies))]})
    return bootdata,dfresults

def resanity(experiments,path='.',**kwargs):
    path=Path(path)
    import numpy.linalg as la
    import seaborn as sns
    tab=pd.DataFrame([])
    figs={}

    ### FROM FINAL MATRIX
    for exp in experiments:
        if 'random' in exp:
            continue
        df, bootdata, results, figdatas,expsimus = experiments[exp]
        alpha_details = hyperplane(df[df.richness>1],list(bootdata['species'].values[0]), etamode=0)[-2]
        Khyp = alpha_details["K"]
        species=list(bootdata['species'].values[0])
        for figdata in figdatas:
            try:
                eta=figdata['eta_sort']
            except:
                etamean=np.mean([np.mean(x,axis=1) for x in results['eta'].values],axis=0)
                compo=list(results['composition'].values[0])
                eta=etamean[ [compo.index(z) for z in figdata['species_sort']]]
            beta,runavg,theomat=np.array(figdata['beta_sort']),np.array(figdata['runavgmat_sort']),np.array(figdata['theomat_sort'])
            for z in (beta,runavg,theomat):
                z[np.isnan(z)]=0
            np.fill_diagonal(beta,0)
            etapred=1-np.dot(beta,eta)
            np.fill_diagonal(beta,1),np.fill_diagonal(runavg,1),np.fill_diagonal(theomat,1)
            ones=np.ones(eta.shape)
            # code_debugger()
            Kmono,Khyp=np.median([list(x) for x in bootdata['Kmono'].values],axis=0),Khyp.values
            Kmono,Khyp=[KK[ [species.index(z) for z in figdata['species_sort']]] for KK in (Kmono,Khyp)]
            dic={
                'Kmono':Kmono,'Khyper':Khyp ,
                'exp':exp.split('_')[0], 'eta':eta,'etapred':etapred,'b.e':np.dot(beta,eta),'r.e':np.dot(runavg,eta),'t.e':np.dot(theomat,eta),
             'nbe':la.norm(np.dot(beta,eta)-ones),'nre':la.norm(np.dot(runavg,eta)-ones),'nte':la.norm(np.dot(theomat,eta)-ones) }
            tab=pd.concat([tab,pd.DataFrame(dic)],ignore_index=1)
    fig, axes = plt.subplots(ncols=2,figsize=np.array(mpfig.rcParams['figure.figsize'])*(1.3,1)  )
    figs['betaeta']=fig
    tab=tab.sort_values('exp')
    tab['exppos']=np.concatenate([[ 1./gp['b.e'].std() for z in range(gp.shape[0])] for exp,gp in tab.groupby('exp')] )
    tab.sort_values('exppos')
    sns.boxplot(x='b.e',y='exp',data=tab,ax=axes[0]),axes[0].set_title('From final matrix'),axes[0].axvline(1)

    if kwargs.get('replicaplots'):
        ### PER REPLICA
        reptab=pd.DataFrame([])
        for exp in experiments:
            if 'random' in exp:
                continue
            df, bootdata, results, figdatas,expsimus = experiments[exp]
            for idx in range(min(results.shape[0],bootdata.shape[0])):
                data,result=bootdata.iloc[idx],results.iloc[idx]
                eta=np.median(result['eta'],axis=1)
                beta,runavg,theomat=np.array(data['beta']),np.array(result['runavgmat']),np.array(result['theomat'])
                for z in (beta,runavg,theomat):
                    z[np.isnan(z)]=0
                cidx=[list(data['species']).index(i) for i in result['composition']]
                beta,theomat,runavg=[c[np.ix_(cidx,cidx)] for c in (beta,theomat,runavg)]
                np.fill_diagonal(beta,0)
                etapred=1-np.dot(beta,eta)
                np.fill_diagonal(beta,1),np.fill_diagonal(runavg,1),np.fill_diagonal(theomat,1)
                ones=np.ones(eta.shape)
                dic={'exp':exp.split('_')[0],'replica':idx, 'eta':eta,'etapred':etapred,'b.e':np.dot(beta,eta),'r.e':np.dot(runavg,eta),'t.e':np.dot(theomat,eta),
                 'nbe':la.norm(np.dot(beta,eta)-ones),'nre':la.norm(np.dot(runavg,eta)-ones),'nte':la.norm(np.dot(theomat,eta)-ones) }
                reptab=pd.concat([reptab,pd.DataFrame(dic)],ignore_index=1)
        plt.suptitle('beta.eta')
        reptab=reptab.sort_values('exp')
        reptab['exppos']=np.concatenate([[ 1./gp['b.e'].std() for z in range(gp.shape[0])] for exp,gp in reptab.groupby('exp')] )
        reptab.sort_values('exppos')
        sns.boxplot(x='b.e',y='exp',data=reptab,ax=axes[1]),axes[1].set_title('From per-bootstrap replicas'),axes[1].axvline(1)


    fig, axes = plt.subplots(ncols=2,figsize=np.array(mpfig.rcParams['figure.figsize'])*(1.3,1)  )
    figs['sanity']=fig
    xs,ys=tab['etapred'], tab['eta']
    rg=np.min(np.concatenate([xs,ys])),np.max(np.concatenate([xs,ys]))
    axes[0].scatter(xs,ys),axes[0].set_xlabel(r"Predicted $\eta$"),axes[0].set_ylabel(r"Observed $\eta$")
    axes[0].plot(rg,rg,linestyle='--')
    xs,ys=tab['Kmono'],tab['Khyper']
    rg=np.min(np.concatenate([xs,ys])),np.max(np.concatenate([xs,ys]))
    axes[1].scatter(xs,ys),plt.xlabel(r"Monoculture $K$"),plt.ylabel(r"Polyculture intercept")
    axes[1].plot(rg,rg,linestyle='--')
    # plt.show()
    plt.suptitle("Lotka-Volterra checks")

    rpath=Path(path+kwargs.get('fname','results'))
    rpath.mkdir()
    if kwargs.get('save',None):
        for fig in figs:
            figs[fig].savefig(rpath+fig+'.svg' )



def get_df(exp):
    df= pd.read_csv(exp+'.csv',sep='\t')
    df['composition']=[ast.literal_eval(c) for c in df.composition]
    return df

def mainboot():
    experiments = {}
    fname = ''
    exps = []
    for s in sys.argv:
        if 'simu' in s:
            exps.append(s)
        if 'data=' in s:
            exps.append(s.split('=')[-1])

    minusbm = 'bm' in sys.argv  ################        MINUSBM

    exps = list(exps)
    for exp in tuple(exps):
        if 'simu' in exp:
            exps.remove(exp)
            nsimus = 1
            if 'random' in exp:
                nsimus = 5
            if 'test' in exp or ('cc' in exp and not 'ccopt' in exp) :
                nsimus = 1
            for n in range(nsimus):
                exps.append('{}_{}'.format(exp, n))
    suffixes={}

    if 'export' in sys.argv:
        export(exps)

    for exp in exps:
        print '=============', exp

        path,suffix=get_path(exp,return_suffix=1)

        suffixes[exp] = suffix
        runbins = ifelse('simu' in exp, 4, 5)
        if not 'simu' in exp:
            df = get_df(exp)
            species = get_species(df)
            if 'detrendK' in sys.argv:
                mono = df[df.richness == 1]
                trendK = {yr: np.mean([cp[c[0]].mean() for c, cp in gp.groupby('composition')]) for yr, gp in
                          mono.groupby('year')}
                df[species] = [c[1:] / trendK[c[0]] for c in df[['year'] + list(species)].values]
            elif  'detrendblock' in sys.argv and 'block' in df:
                if 'total' in df:
                    trend={tuple(yr):np.sum(gp.groupby('plot').mean()['total']) for yr,gp in df.groupby(['year','block'])}
                else:
                    trend = {tuple(yr): np.sum(gp.groupby('plot').mean()[species].values) for yr, gp in
                          df.groupby(['year','block'])}
                df[species] = [c[2:] / trend[(c[0],c[1])] for c in df[['year','block'] + list(species)].values]
            elif 'detrend' in sys.argv:
                if 'total' in df:
                    trend={yr:np.mean(gp.groupby('plot').mean()['total']) for yr,gp in df.groupby('year')}
                else:
                    trend = {yr: np.mean([np.sum(cp[list(c) ].mean()) for c, cp in gp.groupby('composition')]) for yr, gp in
                          df.groupby('year')}
                df[species] = [c[1:] / trend[c[0]] for c in df[['year'] + list(species)].values]
            if 'showtraj' in sys.argv:
                debug_analyze(species, df)

            Nfilter = 'drop2'
            if 'drop5' in suffixes[exp]:
                Nfilter = 'drop5'

            print '  NFILTER:  ', Nfilter

            if 'all' in sys.argv:
                experiments[exp] = (pd.DataFrame({}),)
            else:
                experiments[exp] = (df,)

            if 'bootdata.json' in os.listdir(path) and not 'rerun' in sys.argv:
                bootdata, results = pd.read_json(path + 'bootdata.json'), pd.read_json(path + 'bootresults.json')
                for d in (bootdata, results):
                    for key in d:
                        if key in ('composition', 'species'):
                            d[key] = [tuple(c) for c in d[key]]
                            continue
                        test = d[key].values[0]
                        if hasattr(test, 'keys'):
                            pass
                        elif hasattr(test, '__iter__') and not (
                            isinstance(test, tuple) or isinstance(test, basestring)):
                            # print key, test
                            d[key] = [np.array(m, dtype='float') for m in d[key].values]
            else:
                bootdata, results = bootstrap(df, filter_by=Nfilter, hyperoptions='',
                                              scaling='beta', take_mean_beta=0,
                                              ntrials=40,
                                              testmean='all', infermean='all',
                                              split_mono='monosplit' in sys.argv,
                                              split_criterion= 'nofull',
                                              betalim=ifelse('betalim' in sys.argv,2,15),runbins=runbins)
                if not 'bug' in sys.argv:
                    bootdata.to_json(path + 'bootdata.json'), results.to_json(path + 'bootresults.json')
                    f=open(path+'options','w')
                    f.write(' '.join(sys.argv) )
                    f.close()
        else:
            if 'ccopt' in exp:
                betamode, etamode = 'cc', 'gauss'
            elif 'ccpert' in exp:
                betamode = 'ccpert'
                etamode = 'cc'
            elif 'cc' in exp:
                betamode = etamode = 'cc'
            elif 'pow' in exp:
                betamode, etamode = 'gauss', 'pow'
            else:
                betamode, etamode = 'gauss', 'gauss'
            if 'bootdata.json' in os.listdir(path) and not 'rerun' in sys.argv:
                bootdata, results = pd.read_json(path + 'bootdata.json'), pd.read_json(path + 'bootresults.json')
                for d in (bootdata, results):
                    for key in d:
                        if key in ('composition', 'species'):
                            d[key] = [tuple(c) for c in d[key]]
                            continue
                        test = d[key].values[0]
                        if hasattr(test, 'keys'):
                            pass
                        elif hasattr(test, '__iter__') and not (
                            isinstance(test, tuple) or isinstance(test, basestring)):
                            # print key, test
                            d[key] = [np.array(m, dtype='float') for m in d[key].values]
            else:
                Simu = 8
                if 'S16' in exp:
                    Simu = 16
                elif 'S25' in exp:
                    Simu = 25
                elif 'S50' in exp:
                    Simu = 50
                Simu=ifelse('simu' in etamode,int(Simu*1.3),Simu)
                bootdata, results = simus(betamode=betamode, etamode=etamode, bm=3. / Simu, bs=.66 / np.sqrt(Simu),
                                              ncompo=1, niter=1, S=Simu,runbins=runbins)

                bootdata.to_json(path + 'bootdata.json'), results.to_json(path + 'bootresults.json')
            experiments[exp] = (pd.DataFrame({}),)
            df=None
            species=bootdata['species'].values[0]


        path = Path('data' + ifelse(minusbm, '_bm', '') + '/' + exp + suffixes[exp])
        if not 'composition' in results :#or len(results['composition'].unique().values)<2:
            ls=[(tuple(bootdata['species'].values[0]),results) ]
            testcompos=[None]
        else:
            testcompos = sorted(set([tuple(c) for c in results['composition'].values]))
            figdatas = []
            ls=[(x,y) for x,y in results.groupby('composition')]
        sim=[]
        for comp, gp in ls:
            if len(testcompos) == 1:
                lname = fname
                title = exp
            else:
                print '    COMPO', comp
                suffix = 'comp{}'.format(testcompos.index(comp))
                lname = fname + suffix
                title = exp + suffix

            figdata = None
            noshow = 'noshow' in sys.argv 
            if noshow:
                try:
                    figdata = pickle.load(open(Path(path + lname) + 'figdata.pickle', 'r'))
                except:
                    pass
            if figdata is None:
                # code_debugger()
                figdata = show_bootstrap(bootdata, gp, title=title, runbins=runbins,
                                         path=Path(path),
                                         fname=lname, save=1, debug=0, hold=1, minusbm=minusbm,
                                         redo='rerun' in sys.argv or 'reshow' in sys.argv)
                if noshow:
                    plt.close('all')
            figdatas.append(figdata)
            sim.append(exp_simus(exp,gp,tuple(comp),path=Path(path+lname),runbins=runbins))
        experiments[exp] = experiments[exp] + (bootdata, results, figdatas,sim)


    scorepath=path
    patdat = pattern_score(experiments, path=scorepath,fname=fname, save=not 'simu' in exp,merge='merge' in sys.argv)
    patdat.to_json(Path('data' + ifelse(minusbm, '_bm', '') + '/') + 'pattern_score.json')

    if not 'simu' in exp:
        resanity(experiments, path=Path(path),fname=fname, save=not 'simu' in exp,)
    plt.show()




if __name__ =='__main__':
    import pickle
    for a in sys.argv:
        if 'make=' in a:
            df,ground=make_groundtruth(mu=.35,sigma=.2, replicas=2)
            df.to_csv(a.split('=')[-1]+'.csv',index=None,sep='\t')
    mainboot()
