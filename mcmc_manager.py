import os
import sys, getopt


if __name__ == "__main__":
    index = 0
    filter = 'g'
    npzfolder='/home/dbrout/pysmp_standalone/specnpzfiles'
    outpath = '/home/dbrout/pysmp_standalone/specfitout'
    corioutpath = '/home/dbrout/pysmp_standalone/corispec/'

    try:
        args = sys.argv[1:]
        opt, arg = getopt.getopt(
            args, "i",
            longopts=["index=",'npzfolder=','outpath=','filter='])

    except getopt.GetoptError as err:
        print "No command line arguments"

    for o,a in opt:
        if o in ["-h","--help"]:
            print __doc__
            sys.exit(0)
        elif o in ["-i","--index"]:
            index = a
        elif o in ["--npzfolder"]:
            npzfolder = a
        elif o in ["--outpath"]:
            outpath = a
        elif o in ["--filter"]:
            filter = a

    import scipy.signal
    import numpy as np
    import time

    print 'index',index,

    npzlist = np.asarray(os.listdir(npzfolder),dtype='str')
    newnpzlist = []
    numepochs = []
    for n in npzlist:
        try:
            if not os.path.exists(corioutpath + '/' + n.split('.')[0] + '.smp'):
                #if not os.path.exists(corioutpath + '/' + n.split('.')[0] + '.mcmcout'):
                numepochs.append(np.load(npzfolder+'/'+n)['Nimage'])
                newnpzlist.append(n)
        except:
            pass
    numepochs = np.asarray(numepochs)
    npzlist = np.asarray(newnpzlist,dtype='str')

    npzlist = npzlist[np.argsort(numepochs)]
    npzfile = npzlist[int(index)]
    #npzfile = os.listdir(npzfolder)[int(index)]
    inp = np.load(npzfolder+'/'+npzfile)
    print inp.keys()
    print 'numepochs', inp['Nimage']

    #outpath = 'fitout/'
    lcout = outpath+'/'+npzfile.split('.')[0]
    chainsnpz = outpath+'/'+npzfile.split('.')[0] + '_chains.npz'
    stdoutfile = outpath+'/'+npzfile.split('.')[0] + '.log'
    smpfile = outpath+'/'+npzfile.split('.')[0] + '.smp'
    if os.path.exists(corioutpath+'/'+npzfile.split('.')[0] + '.smp'):
        print 'already ran'
    else:
        import mcmc

        #print inp['mjdflag']
        #raw_input()
        a = mcmc.metropolis_hastings(
            galmodel=inp['galmodel']
            , modelvec=inp['modelvec']
            , galstd=inp['galstd']
            , modelstd=inp['modelstd']
            , data=inp['data']
            , psfs=inp['psfs']
            , weights=inp['weights']
            , substamp=inp['substamp']
            , Nimage=inp['Nimage']
            , maxiter = 1500000 #, maxiter=inp['maxiter']
            , mask=inp['mask']
            , sky=inp['sky']
            , mjd=inp['mjd']
            , gewekenum= 5000 #, gewekenum=inp['gewekenum']
            , skyerr=inp['aper_skyerr']
            , useskyerr=inp['useskyerr']
            , usesimerr=inp['usesimerr']
            , flags=inp['smpdictflag']
            , fitflags=inp['fitflags']*0.
            , psf_shift_std=inp['psf_shift_std']
            , xoff=0.
            , yoff=0.  # .06
            , shiftpsf=inp['shiftpsf']
            , fileappend=inp['fileappend']
            , stop=False
            , skyerr_radius=inp['skyerr_radius']
            , outpath=outpath
            , compressionfactor=inp['compressionfactor']
            , fix_gal_model=inp['fix_gal_model']
            , pixelate_model=None
            , burnin=inp['burnin']
            , lcout=lcout
            , chainsnpz=chainsnpz
            , mjdoff=inp['mjdoff']
            , dontsavegalaxy=True
            , log=None
            , isfermigrid=False
            , isworker=False
            , x=inp['x']
            , y=inp['y']
            , psffile=inp['psffile']
            #, psfcenter=inp['psfcenter']
            , model_errors=True
            , survey=inp['survey']
            , fileroots=inp['fileroots']
            , scalefactor=inp['scalefactor']
            , gain=inp['gain']
            , dobkg=False
            , sigmazpt=inp['sigmazpt']
            , fakemag=inp['fakemag']
            , fitzpt=inp['fitzpt']
            , fitzpterr=inp['fitzpterr']
            , fakezpt=inp['fakezpt']
            , datafilenames=inp['datafilenames']
            , nightlyoffx=inp['nightlyoffx']
            , nightlyoffy=inp['nightlyoffy']
            , sstime=inp['sstime']
            , stdoutfile = None
            , smpfile = smpfile
            , peakmjd=inp['peakmjd']
            , idobs=inp['idobs']
            , idcoadd=inp['idcoadd']
            , diffim_flux=inp['diffim_flux']
            , diffim_fluxerr=inp['diffim_fluxerr']
            , ra=inp['ra']
            , dec=inp['dec']
            , smpdictflag=inp['smpdictflag']
            , mjdflag=inp['mjdflag']
            , descriptiveflag=inp['descriptiveflag']
            , rmsaddin=inp['rmsaddin']
            , gewekediag=inp['gewekediag']
            , imfilename=inp['imfilename']
            , weightfilename=inp['weightfilename']
            , zptfilename=inp['zptfilename']
            , filt = inp['filt']

        )

    # modelvec, modelvec_uncertainty, galmodel_params, \
    # galmodel_uncertainty, modelvec_nphistory, galmodel_nphistory, \
    # sims, xhistory, yhistory, accepted_history, pix_stamp, \
    # chisqhist, redchisqhist, stamps, chisqs, ndof, gewekediag, \
    # covar, corr = a.get_params()












