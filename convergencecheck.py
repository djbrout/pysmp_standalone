import numpy as np
import scipy.signal
import scipy.misc
import scipy.ndimage as nd
# from . import Matplot, flib
# from .utils import autocorr, autocov
from copy import copy
from numpy import corrcoef, sum, log, arange
from numpy.random import rand
# from pylab import pcolor, show, colorbar, xticks, yticks
# import pylab as plt
import time
import os
import math
import matplotlib as m

m.use('Agg')
import matplotlib.pyplot as plt
import scipy.interpolate as interpol
import scipy.ndimage
import dilltools as dt
from matplotlib.backends.backend_pdf import PdfPages
import gc
# import chkpsf_fast
import matplotlib.mlab as mlab
import math
import build_psfex
import scipy.stats
import geweke as g
# import pymc
# import rdpsf
import sys
from scipy.fftpack import fft, ifft, fft2, ifft2
import multiprocessing
import time
import os
# import pyfftw
from numpy import corrcoef, sum, log, arange
from numpy.random import rand
from pylab import pcolor, show, colorbar, xticks, yticks


def check_geweke(chain,burnin=.3):

    num_iter = int(round(len(chain) * (1. - burnin)))
    start_iter = int(round(len(chain) * (burnin)))
    if num_iter < 200:
        print 'num iter too small'
        return True
    hasnotconv = False

    #self.gewekediag = np.zeros_like(self.modelstd)

    #self.modelvec_nphistory = np.zeros((num_iter, len(self.modelvec)))
    # print num_iter, start_iter
    # print len(self.modelvechistory)
    #for i in np.arange(num_iter - 1):
    #    self.modelvec_nphistory[i, :] = self.modelvechistory[int(i + start_iter)]

    # print num_iter
    # print start_iter
    # print len(self.modelvec_nphistory[0, :])
    #for param in range(len(self.modelvec_nphistory[0, :])):
        # # print self.modelvec_nphistory.shape
        # # print self.modelvec_nphistory[param,:].shape
        # # print len(np.unique(self.modelvec_nphistory[param,:]))
        # # print len(np.unique(self.modelvec_nphistory[:,param]))
        # # if len(np.unique(self.modelvec_nphistory[:, param])) == 1:
        # #     print 'asdf'
        # #     self.gewekediag[param] = -999.
        # #     continue
        # if self.modelstd[param] == 0:
        #     # print 'stdzero'
        #     self.gewekediag[param] = -999.
        #     continue

        #try:
    gw = g.geweke(chain, intervals=2, first=.4, last=.5)
    print gw
    raw_input()
    # gew = []
    # for gg in gw:
    #     gew.append(gg[1])
    # gew = np.array(gew)
    #except:
    # #        gew = np.array([999., 999.])
    #
    #     # print gew.shape
    #     # raw_input('gews')
    #     # gew = dt.geweke(self.modelvec_nphistory[:, param])
    #
    #  #   self.gewekediag[param] = np.mean(gew)
    #
    #     # self.gewekediag[param] = np.mean(np.abs(geweke[:, 1]))
    #
    #     print param, gew
    #
    #     if np.any(np.abs(gew) > 2.):
    #         msg = "Epoch %s has not properly converged" % param
    #         # if assert_:
    #         #     raise AssertionError(msg)
    #         # else:
    #         print(msg)
    #
    #         hasnotconv = True
    #         # return True

    return gw











def getgeweke(chain,burnin=.3):
    start_iter = int(round(len(chain) * (burnin)))
    try:
        gw = g.geweke(chain[start_iter:], intervals=1, first=.4, last=.5)
        # gew = []
        # for gg in gw:
        #     gew.append(gg[1])
        # gew = float(gew)
    except:
        gw = np.nan
    return gw



outdir = '/project/projectdirs/dessn/dbrout/simv3.0/convergencetest/'

allchains = []
for f in os.listdir(outdir)[:40]:
    if '_chains.npz' in f:
        print f
        allchains.append(np.array(np.load(outdir+f)['local_galchain']))

plt.figure(figsize=(12,9))
for j,chain in enumerate(allchains):
    gvec = []

    for i in range(len(chain)):
        print i,getgeweke(chain[:i])
        raw_input()
        gvec.append(getgeweke(chain[:i]))

    print j,'of',len(allchains)
    plt.plot(np.arange(len(chain))*100.,gvec

plt.savefig('convergence.png')


















