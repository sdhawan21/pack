"""
Gaussian Process light curve fits for SNIa (NIR, t2 fits )

To use on CSP (Contreras+2010 & Stritzinger+2011), but also for the upcoming CfA data release (Friedman+ 2014)

Based on george (Dan Foreman-Mackey)
"""

#package that contains the proper light curve reading functions

from pack import dist	


import emcee 
import numpy as np
import matplotlib.pyplot as plt
import sys
import george

from george.kernels import ExpSquaredKernel

insn=sys.argv[1]
inband=sys.argv[2]
pt='/Users/lapguest/all_paper/files_snpy/'
def set_arr(insn, inband):
	"""
	Read the light curve 
	"""	
	lc=dist.rd_lc(insn, inband)	
	return lc
#define model, likelihood, prior and posterior if you want to do a full Bayesian treatment 

#good for correlated noise. If you don't have correlated noise, do yourself a favour and skip to the main function and do the simple model fit (sorry if that sounds rude, really sorry, ill buy you a beer. hey, nobody else is going to read it, whats the point, be rude to yourself)

def model1(params, t):
    m, b, amp, loc, sig2 = params
    return m*t + b + amp * np.exp(-0.5 * (t - loc) ** 2 / sig2)

def lnlike1(p, t, y, yerr):
    return -0.5 * np.sum(((y - model1(p, t))/yerr) ** 2)

def lnprior1(p):
    m, b, amp, loc, sig2 = p
    if (-10 < m < 10 and  -10 < b < 10 and -10 < amp < 10 and
            -5 < loc < 5 and 0 < sig2 < 3):
        return 0.0
    return -np.inf

def lnprob1(p, x, y, yerr):
    lp = lnprior1(p)
    return lp + lnlike1(p, x, y, yerr) 




def main():
	"""
	Meat of the code

	"""

	par=np.loadtxt(pt+'tmax_dm15.dat', dtype='string')

	tbmax=float(par[par[:,0]==insn][0][1])
	
	lc=set_arr(insn, inband)
	
	ph=lc['MJD']-tbmax
	#condition for second maximum 
	##TODO: GUI for selecting region
	cond=(ph>=10.0) & (ph<=40.0)
	#define the data in the region of interest

	ph1=ph[cond]
	
	mag=lc[inband][cond]

	magerr=lc['e_'+inband][cond]
	
	# Set up the Gaussian process.
	kernel = ExpSquaredKernel(1e2)
	gp = george.GP(kernel)
	
	# Pre-compute the factorization of the matrix.
	gp.compute(ph1, magerr)

	#vertically stack it for no apparent reason 
	arr=np.vstack([ph1, mag, magerr]).T
	"""
	initial = np.array([0, 0, -1.0, 0.1, 0.4])
	ndim = len(initial)
	nwalkers=100

	p0 = [np.array(initial) + 1e-8 * np.random.randn(ndim)
      	for i in xrange(nwalkers)]
	sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob1, args=(ph1, mag ,magerr))

	print("Running burn-in...")
	p0, _, _ = sampler.run_mcmc(p0, 500)
	sampler.reset()

	print("Running production...")
	sampler.run_mcmc(p0, 1000)
	"""
	#define x array for applying the fit
	t = np.linspace(ph1.min(), ph1.max(), 500)

	mu, cov = gp.predict(mag, t)
	
	#condition for peak
	mpeak=(mu==min(mu))
	
	#calculate the standard deviation from covariance matrix 
	std = np.sqrt(np.diag(cov))
	
	print t[mpeak][0], min(mu), std[mpeak][0]

	plt.errorbar(ph1, mag, magerr, fmt=".k", capsize=2, label='data')
	
	#plt.plot(t, mu, 'k:', label='best fit')
	plt.fill_between(t, mu-std, mu+std)
	plt.legend(loc=0)	
	plt.ylim(plt.ylim()[::-1])
	plt.show()
if len(sys.argv)==3:
	main()
else:
	print "Usage: python"+sys.argv[0]+"<SN> <band> "

