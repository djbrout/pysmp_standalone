import os
import sys, getopt
import time

if __name__ == "__main__":
    index = 0
    filter = 'g'
    npzfolder='/home/dbrout/pysmp_standalone/specnpzfiles'
    outpath = '/home/dbrout/pysmp_standalone/specfitout'
    corioutpath = '/home/dbrout/pysmp_standalone/corispec/'
    npzfile = None

    try:
        args = sys.argv[1:]
        opt, arg = getopt.getopt(
            args, "i",
            longopts=["index=",'npzfolder=','outpath=','filter=','sn='])

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
        elif o in ['--sn']:
            npzfile = a

    import scipy.signal
    import numpy as np
    import time

    print 'index',index,

    #time.sleep(float(index))

    npzlist = np.asarray(sorted(os.listdir(npzfolder)),dtype='str')
    newnpzlist = []
    numepochs = []
    for n in npzlist:
        try:
            if not os.path.exists(outpath + '/' + n.split('.')[0] + '.smp'):
                #if not os.path.exists(outpath + '/' + n.split('.')[0] + '.running'):
                if np.load(npzfolder+'/'+n)['peakmjd'] > 0:
                    numepochs.append(np.load(npzfolder+'/'+n)['Nimage'])
                    newnpzlist.append(n)
        except:
            pass
    numepochs = np.asarray(numepochs)
    npzlist = np.asarray(newnpzlist,dtype='str')

    #npzlist = npzlist[np.argsort(numepochs)]
    #npzfile = npzlist[int(index)]
    if npzfile is None:
        npzfile = npzlist[int(index)]

    #smprunningfile = outpath+'/'+npzfile.split('.')[0] + '.running'
    #os.system('touch '+smprunningfile)

    #npzfile = os.listdir(npzfolder)[int(index)]
    inp = np.load(npzfolder+'/'+npzfile)

    #outpath = 'fitout/'
    lcout = outpath+'/'+npzfile.split('.')[0]
    chainsnpz = outpath+'/'+npzfile.split('.')[0] + '_chains.npz'
    stdoutfile = outpath+'/'+npzfile.split('.')[0] + '.log'
    smpfile = outpath+'/'+npzfile.split('.')[0] + '.smp'

    if os.path.exists(outpath+'/'+npzfile.split('.')[0] + '.smp'):
        print 'already ran'
    else:
        import mcmc
        print npzfile
        print inp['compressionfactor']
        #raw_input()

        #raw_input('asdf')
        a = mcmc.metropolis_hastings(
            galmodel=inp['galmodel']
            , modelvec=inp['modelvec']
            , galstd=inp['galstd']/1.1
            , modelstd=inp['modelstd']/1.1
            , data=inp['data']
            , psfs=inp['psfs']
            , weights=inp['weights']
            , substamp=inp['substamp']
            , Nimage=inp['Nimage']
            , maxiter = 2000000 #, maxiter=inp['maxiter']
            , mask=inp['mask']
            , sky=inp['sky']
            , mjd=inp['mjd']
            , gewekenum= 5000 #, gewekenum=inp['gewekenum']
            , skyerr=inp['aperskyerr']
            , useskyerr=inp['useskyerr']
            , usesimerr=inp['usesimerr']
            , flags=inp['smpdictflag']
            , fitflags=inp['fitflags']*0.
            , psf_shift_std=0.0002
            , xoff=0.
            , yoff=0.  # .06
            , shiftpsf=True
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
            #, mjdoff=inp['mjdoff']
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
            , stdoutfile = stdoutfile
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












