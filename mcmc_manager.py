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

    r = np.random.random((64,64))
    t1 = time.time()
    for i in range(10000):
        a = scipy.signal.fftconvolve(r,r*5.)
    t2  = time.time()
    print 'index',index,'time',t2-t1
