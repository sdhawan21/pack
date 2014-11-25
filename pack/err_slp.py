import numpy as np
import matplotlib.pyplot as plt
import sys
import statsmodels.api as sm


from scipy.interpolate import interp1d
from scipy.odr import *
from scipy.stats import pearsonr
from sklearn import linear_model

#define path
pt='/home/sdhawan/tests_paper/csp_sn/sec_max_files/'


#input variables
#band, t_input, error_t
band=sys.argv[1]
t=float(sys.argv[2])
et=float(sys.argv[3])

#print band

#input files 
ufl=np.loadtxt('../tables/u_flags.txt', dtype='string', skiprows=1)
if sys.argv[5]=='fl':
	ufl=ufl[ufl[:,-1]=='Y']
tj=np.loadtxt(pt+band+'_sec_max_csp.dat', dtype='string')
y2j=np.loadtxt('../y2j.txt', dtype='string')
y2j=np.loadtxt(pt+'y_sec_max_csp.dat', dtype='string')

#shorthand for function
rn=np.random.normal
fac=float(sys.argv[6])
def spl_fit(arr, val):
	"""
	interpolate value for Nickel mass
	"""
	real=rn(arr[:,0], arr[:,1])
	terp=interp1d(real, rn(arr[:,2], arr[:,3]), kind='cubic')
	l=np.linspace(real.min(), real.max())
	gpl=terp(l)
	l1=abs(l-val)
	return gpl[l1==min(l1)][0]
	
def dist_fac(h0):
	return (70/float(h0))**2
def arr_crt(ufl, tj):
	lb2=np.array([float(i[1]) for i in ufl if i[0] in tj[:,0]])
	elb2=np.array([float(i[2]) for i in ufl if i[0] in tj[:,0]]) 
	jt2=np.array([float(tj[tj[:,0]==i[0]][0][1]) for i in ufl if i[0] in tj[:,0]])
	jet2=np.array([float(tj[tj[:,0]==i[0]][0][2]) for i in ufl if i[0] in tj[:,0]])
	return lb2, elb2, jt2, jet2
def lbol_red(tval):
	"""
	function to estimate Lbol from t_2
	"""
	
	#to evaluate the effect of hubble constant values for the distance measurement of the calibration sample
	ufl1=ufl[ufl[:,-1]=='Y']
	ufl2=ufl[ufl[:,-1]=='N']
	
	#print ufl1[:,0]
	#define arrays for x, y, ex, ey (not in that order)
	lb2, elb2, jt2, jet2=arr_crt(ufl1, tj)
	
	
	lb2*=dist_fac(fac)
	
	lb3, elb3, jt3, jet3=arr_crt(ufl2, tj)
	#array for the y2j conversion (if specified by command line)
	if int(sys.argv[4])==1:
		lb1, elb1, iyt2, iyet2=arr_crt(ufl1, y2j)
		lb4, elb4, ijt2, ijet2=arr_crt(ufl2, y2j)
		lb1*=dist_fac(fac)	
		lb=np.concatenate([lb1, lb2, lb3, lb4]); elb=np.concatenate([elb1, elb2, elb3, elb4]); t2=np.concatenate([ iyt2, jt2, jt3, ijt2])
		et2=np.concatenate([ iyet2, jet2, jet3, ijet2])
	
	else:
		lb=np.concatenate([lb3, lb2]); elb=np.concatenate([elb3, elb2]); t2=np.concatenate([jt3, jt2])
		et2=np.concatenate([jet3, jet2])
		#lb=lb2; elb=elb2; t2=jt2; et2=jet2
	
	
		
	#print len(lb)
	#load constructor for fit
	clf=linear_model.LinearRegression()
	
	#x=[jt2, yt2]
	#x=sm.add_constant(x)
	#y=np.array([lb,lb1])
	#x=np.array([jt2, yt2])
	#clf.fit(x, y)
	
	#print clf.coef_
	#return 0
	
	#define Nsamp
	nsamp=100
	
	
	
	#perform least-squares
	rd=RealData(t2, lb, sx=et2, sy=elb)
	def f(B,x):
		return B[0]*x+B[1]
	f=Model(f)
	out=ODR(rd, f, beta0=[1., 2.])
	o=out.run()
	val=o.beta
	err=o.sd_beta
	#print "The total number of NS, errs, best fit values and pearsonr values are"
	#print "The slope and intercept with their errors are: \t", val, err, pearsonr(t2, lb)
	
	
	arr=np.vstack([t2, et2, lb, elb]).T
	#np.savetxt('../out_files/bivar_regress.txt', arr, fmt='%s')
	""" 
	spline fits 
	"""
	ests=[]
	for k in range(nsamp):
		try:
			ests.append(spl_fit(arr, tval[0]))
		except:
			tval
			
			
	#do monte carlo for Mni (using Arnett+fixed rise)
	ar=[(np.random.normal(val[0], err[0])*np.random.normal(tval[0], tval[1])+np.random.normal(val[1], err[1])) for i in range(3000)]

	mx=val[0]*tval[0]; c=val[1]
	return ar#mx+c, mx*((tval[1]/tval[0])+(err[0]/val[0])+err[1]), #np.array(ar)


def coeff(dm15):
	"""
	calculate coefficient of arnett's relation
	"""
	rt=17.5-5*(dm15-1.1)
	lm=6.45e43*np.exp(-rt/8.8)+1.45e43*np.exp(-rt/111.1)
	return lm/1e43
def main():
	v=t
	ev=et                              
	ar=lbol_red([v, ev])
	print np.mean(ar), np.std(ar), fac#, min(ar), max(ar)#np.mean(ar)/2, (np.mean(ar)/2)*((np.std(ar)/np.mean(ar))+(0.3/2.0))
	
if len(sys.argv)==7:
	main()
else:
	print "Usage: python"+sys.argv[0]+"<band> <t2 (obj)> <et2(obj)> <make combine> <flags> <factors on errors>"
