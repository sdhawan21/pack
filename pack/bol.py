from scipy.interpolate import interp1d

import numpy as np

rn=np.random.normal
class bol_func:
	"""
	Class to store all functions required for evaluating parameters from bolometric light curves (calculated using mklcbol.pro)

	Functions also to calculate the rise from half maximum a la contardo 2000


	"""

	def bolpeak(self, fil):
		"""
		Interpolate to find bolometric peak ##TOADD: exceptions for when the peak is not measured
		
		note: written on top of interp1d from scipy.interpolate (boss function it is!)
		"""
		lc=np.loadtxt(fil)
		ter=interp1d(lc[:,0], lc[:,1], kind='cubic')
		l=np.linspace(lc[:,0].min(), lc[:,0].max())
		gpl=ter(l)
		return max(gpl), l[gpl==max(gpl)][0]
	
	def err_peak(self, arr):
		real=rn(arr[:,1], arr[:,2])
		ph=arr[:,0]	
		
		spl=interp1d(ph, real, kind='cubic')
		l=np.linspace(ph.min(), ph.max())
		
		gpl=spl(l)	
		return min(gpl), l[gpl==min(gpl)][0]


	def dm15_bol(self,fil):

		"""

		Delta m 15 for the bolometric light curve 
		"""
		lc1=np.loadtxt(fil)
		
		sp=interp1d(lc1[:,0], lc1[:,1], kind='cubic')
		
		l=np.linspace(min(lc1[:,0]), max(lc1[:,0]), 100)
		
		gp=np.log10(sp(l))
		
		tm=l[gp==max(gp)]
		
		ph=lc1[:,0]-tm
		
		if min(ph)<0 and max(ph)>15:
		
			ph1=abs(l-15)
			m15=gp[ph1==min(ph1)][0]
		
			return max(gp)-m15#np.log10(max(gp))-np.log10(m15)	
		else:
			return 99.0
	def late_decl(self, arr1, ran=[200, 400]):
		arr=arr1[(arr1[:,0]>ran[0]) & (arr1[:,0]<ran[1])]
		if len(arr)>=4:
			A=np.vstack([arr[:,0], np.ones(len(arr))]).T
			m=np.linalg.lstsq(A, arr[:,1])[0]
			return m
		else:
			return [99.0, 99.0]
class mc:
	"""
	Class of MC functions for testing distributions and non-gaussian behaviour
	"""
	def ar_crt(self, n,	t0=[1., 2.]):
		ar=[rn(t0[0], t0[1]) for i in range(n)]
		return ar
	
	def rel_ni(self, mb,n):
		ar=[pow(10, -0.4*(rn(mb[0], mb[1])+rn(19.841, 0.020))) for k in range(n)]
		return ar




class reden:
	
	b14=np.loadtxt("/home/sdhawan/bol_ni_ej/burns14_ebv.tex", dtype='string', delimiter='&')
	def b14_av(self, SN):
		bf=self.b14
		row=bf[bf[:,0]==SN+' '][0]
		try:
			
			ebv=float(row[5][2:7])
			
			rv=float(row[7][2:5])
			
			return ebv*rv
		except:
			return "There is no R_v with this method"
	def b14_sbv(self, SN):
		bf=self.b14
		
		try:
			row=bf[bf[:,0]==SN+' '][0]
			ebv=float(row[3][1:6])
			
				
			return ebv
		except:
			return 99.0	#what a failed attempt at trying to amuse yourself

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
def arn_coef(rt, alp=1):
	"""
	For a given rise time, calculate the coefficient for the relation between Nickel mass and peak bolometric luminosity (arnett's rule, instantaneous energy deposition is output energy at max )
	

	Default alpha is 1 (arguments in Branch+ 1992, stritzinger 2006)
	"""

	eni=6.45e43*np.exp(-rt/8.8)+1.45e43*np.exp(-rt/111.1)

	return alp*eni/1e43

def arn_coef_mc(n, rt=[19, 3], alp=1):
	"""
	n realisations of the coefficient from arnett's rule

	rise time value, default from Stritzinger+

	"""

	ar=[arn_coef(rn(rt[0], rt[1]), alp=1) for k in range(n)]	

	return np.mean(ar), np.std(ar)

def t_half_rise(fil, kind):
	"""
	Calculate the time to rise to max from half the luminosity 
	
	kind: spline, polyfit
	"""
	
	peak=bol_func().bolpeak(fil); tmax=peak[1]
	#load LC
	lc=np.loadtxt(fil)
	#use only the pre-max light curve 
	ph=lc[:,0]-tmax
	ph1=ph[ph<0]; mag1=lc[:,1][ph<0]
	
	if kind == "spline":
		
	
		spl=interp1d(ph1, mag1, kind='cubic')
		l=np.linspace(ph1.min(), ph1.max()); gpl=spl(l)
		arr=abs(gpl-(peak[0]/2.0))
		minval=l[arr==min(arr)][0]
		
		if minval>min(ph):
			return minval
		else:
			return 99.0	

	if kind=="polyfit":
		coef=np.polyfit(ph1, mag1, 2.0)
		lp=np.linspace(-20, 0.0)
		magval=coef[0]*(lp**2)+coef[1]*lp+coef[2]
		mhalf=peak[0]/2.0
		thalf=lp[abs(magval-mhalf)==min(abs(magval-mhalf))]
		
		return thalf
def scal14_ni(x1, n):
	ar=np.array([rn(x1[0], x1[1])*rn(0.100, 0.020)+rn(0.478, 0.023) for k in range(n)])
	return np.mean(ar), np.std(ar)	








