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
    inputs = np.load(npzfile)

    import mcmc

    a = mcmc.metropolis_hastings(*inputs)

    # modelvec, modelvec_uncertainty, galmodel_params, \
    # galmodel_uncertainty, modelvec_nphistory, galmodel_nphistory, \
    # sims, xhistory, yhistory, accepted_history, pix_stamp, \
    # chisqhist, redchisqhist, stamps, chisqs, ndof, gewekediag, \
    # covar, corr = a.get_params()












