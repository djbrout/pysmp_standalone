import os
import sys, getopt


if __name__ == "__main__":
    index = 0
    try:
        args = sys.argv[1:]
        opt, arg = getopt.getopt(
            args, "i",
            longopts=["index="])

    except getopt.GetoptError as err:
        print "No command line arguments"

    for o,a in opt:
        if o in ["-h","--help"]:
            print __doc__
            sys.exit(0)
        elif o in ["-i","--index"]:
            index = a

    


    import scipy.signal
    import numpy as np
    import time

    print 'index',index,

    npzfile = os.listdir('npzfiles/')[int(index)]
    inp = np.load('npzfiles/'+npzfile)

    outpath = 'fitout/'
    lcout = outpath+'/'+npzfile.split('.')[0]
    chainsnpz = outpath+'/'+npzfile.split('.')[0] + '_chains.npz'
    stdoutfile = outpath+'/'+npzfile.split('.')[0] + '.log'
    import mcmc

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
        , maxiter=inp['maxiter']
        , mask=inp['mask']
        , sky=inp['sky']
        , mjd=inp['mjd']
        , gewekenum=inp['gewekenum']
        , skyerr=inp['skyerr']
        , useskyerr=inp['useskyerr']
        , usesimerr=inp['usesimerr']
        , flags=inp['flags']
        , fitflags=inp['fitflags']
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
        , pixelate_model=inp['pixelate_model']
        , burnin=inp['burnin']
        , lcout=lcout
        , chainsnpz=chainsnpz
        , mjdoff=inp['mjdoff']
        , dontsavegalaxy=True
        , log=inp['log']
        , isfermigrid=False
        , isworker=False
        , x=inp['x']
        , y=inp['y']
        , psffile=inp['psffile']
        , psfcenter=inp['psfcenter']
        , model_errors=True
        , survey=inp['survey']
        , fileroots=inp['fileroots']
        , scalefactor=inp['scalefactor']
        , gain=inp['gain']
        , dobkg=False
        , sigmazpt=inp['sigmazpt']
        , fakemag=inp['fakemag']
        , fitzpt=inp['fitzpt']
        , fakezpt=inp['fakezpt']
        , datafilenames=inp['datafilenames']
        , nightlyoffx=inp['nightlyoffx']
        , nightlyoffy=inp['nightlyoffy']
        , sstime=inp['sstime']
        , stdoutfile=stdoutfile

    )

    # modelvec, modelvec_uncertainty, galmodel_params, \
    # galmodel_uncertainty, modelvec_nphistory, galmodel_nphistory, \
    # sims, xhistory, yhistory, accepted_history, pix_stamp, \
    # chisqhist, redchisqhist, stamps, chisqs, ndof, gewekediag, \
    # covar, corr = a.get_params()












