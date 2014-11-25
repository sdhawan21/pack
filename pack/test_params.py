import numpy as np
import sys
import matplotlib.pyplot as plt
import statsmodels.api as sm

from scipy.interpolate import interp1d
from glob import glob
from scipy.stats import pearsonr

#shorthand for function
rn=np.random.normal
#path to dm15 values 
pt="/home/sdhawan/tests_paper/ni/files_snpy/"
pt1="/home/sdhawan/tests_paper/csp_sn/sec_max_files/"
def bolpeak(fil):
	lc=np.loadtxt(fil)
	ter=interp1d(lc[:,0], lc[:,1], kind='cubic')
	l=np.linspace(lc[:,0].min(), lc[:,0].max())
	gpl=ter(l)
	return max(gpl), l[gpl==max(gpl)][0]

def av_err(fil):
	fin=open(fil, 'r')
	ls=[]
	for row in fin:
		ls.append(row.split())
	return [ls[3][1], ls[4][3], ls[14][2]]

def err_bolpeak(fil):
	maxar=[]
	
	#take 1000 realisations of the light curve
	for k in range(1000):
		lc=np.loadtxt(fil)
		real=rn(lc[:,1], lc[:,2]) 
		
		#interpolation routine written within scipy, as Marco says "its fucknig amazing"
	
		ter=interp1d(lc[:,0], real, kind='cubic')
		l=np.linspace(lc[:,0].min(), lc[:,0].max())
		gpl=ter(l)
		maxar.append(max(gpl))
	
	return np.array(maxar)

def ext_name(fil):
	fout=open(fil, 'r')
	ls=[]
	for row in fout:
		ls.append(row.split())
	return ls[3][1]
def boots(arr):
	barr=[]
	for k in range(len(arr)):
		ind=int(np.random.uniform(0, 1)*len(arr))
		barr.append(arr[ind])
	return np.array(barr)
def ni_ddc(val):
	ddc=np.loadtxt('../comp_ddc/lpeak_m56ni.dat', usecols=(7, 1))
	tarr=[]
	for l in range(1000):
		ddc1=ddc
		#ddc1=boots(ddc)
		#ddc1=ddc1[ddc1[:,0].argsort()]
		#ddc1=np.unique(ddc1)
		#print ddc1
		err=np.random.uniform(0.0, 0.1, len(ddc1))
		real_X=rn(ddc1[:,0], np.ones(len(ddc1))*1e-10)
		spl=interp1d(real_X, rn(ddc1[:,1], err) , kind='cubic')
		l=np.linspace(real_X.min(), real_X.max())
		gpl=spl(l)
		l1=abs(l-val)
		tarr.append(gpl[l1==min(l1)][0])
	return np.array(tarr)
def rise_err(dm, lbol):
	"""
	
	calculate error in rise time by monte carlo 
	
	"""
	fac_arr=[]
	
	for k in range(10000):
	
		# Scalzo 2014 equation to cover the dm15-t,Rb locus of Ganeshalingham '11
		t=17.5-5*(rn(dm[0], dm[1])-1.1)
		
		# arnett's rule. Conservative 2 day uncertainty in rise time as in Scalzo 2014
		
		et=float(sys.argv[1])
		
		lb=6.45e43*np.exp(-rn(t, et)/8.8)+1.45e43*np.exp(-rn(t, et)/111.1)
		
		#normalize the coefficient
		fac=lb/1e43
		fac_arr.append(fac)
		
	return np.mean(fac_arr), np.std(fac_arr)

def err(x, y, sx):
	return sx/y 
	
	
def mbol(lbol):
	"""
	Convert bol luminosity to bolometric magnitude using solar references
	
	For comparison with Li+2011 (LOSS sample) and the models of ruiter+2013 on violent mergers and the relation 	   between primary mass and Mbol 
	"""
	return 4.8-np.log10(lbol*1e43/3.84e33)


def nierr_comp():
	"""
	Compare the errors from the calculations from stephan's code to Stritzinger+ (and scalzo+)
	
	"""
	s06=np.loadtxt('../s06_recon/s06_ej_ni.txt', usecols=( 8, 9), delimiter='&', dtype='string')
	s14=np.loadtxt('../s14_files/scal_14_ni.tex', delimiter='&', dtype='string', usecols=(1, 2))
	
	ni=np.array([[float(i[0][1:5]), float(i[0][7:10])] for i in s06])
	scal_ni=np.array([[float(i[1][2:6]), float(i[1][11:15])] for i in s14])
	
	
	eni=lbol[:,0]*((0.3/2.0)+(lbol[:,1]/lbol[:,0]))/2.0
	
	plt.hist(ni[:,1], alpha=0.3, label="S06")
	plt.hist(eni, alpha=0.3, label="evaluated")
	plt.hist(scal_ni[:,1], alpha=0.3, label='Sc14')
	plt.legend(loc=0)
	plt.show()	

def ir_contrib(sn):
	sset=sorted(glob('../lcbol_distrib/'+sn+'*Y*J*.dat'))
	max_arr=[]
	for k in sset:
		inp=np.loadtxt(k)
		plt.plot(inp[:,0], np.log10(inp[:,1]))
		max_arr.append(bolpeak(k)[0])
	max_arr=np.array(max_arr); print max_arr, max_arr.min()/max_arr.max()
	
	#plt.show()

lbol=np.loadtxt('../tables/u_flags.txt', skiprows=1, usecols=(1, 2))
def main():
	
	

	if len(sys.argv)==2:
		inf=sys.argv[1]

	sset=sorted(glob('../lcbol_distrib/finfiles/*.dat'))
	
	t2inp=np.loadtxt(pt1+'j_sec_max_csp.dat', dtype='string')
	
	bp_arr=[]
	
	nmarr=[]; t2arr=[]
	
	sn=sys.argv[1]
	ir_contrib(sn)
	
	
	
	
	return 0
	"""
	for k in sset:
	
		bp=bolpeak(k)
		if ext_name(k) in t2inp[:,0]:
			bp_arr.append(bp[0])
	
			nmarr.append([ext_name(k), bp[0], np.std(err_bolpeak(k)/1e43), float(t2inp[t2inp[:,0]==ext_name(k)][0][1])])
			print ext_name(k)
			t2arr.append(float(t2inp[t2inp[:,0]==ext_name(k)][0][1]))
	"""
	in_arr=np.loadtxt('../out_files/bol_ni.dat', usecols=(1, 2, 3))
	#bp_arr=np.array(bp_arr);t2arr=np.array(t2arr)
	#np.savetxt("../out_files/bol_ni.dat",nmarr, fmt='%s')
	#plt.hist(bp_arr/1e43, alpha=0.3)
	#plt.hist(lbol[:,0], alpha=0.3)
	#plt.show()
	
	
	res=sm.OLS(in_arr[:,0]/2e43, in_arr[:,2]).fit()
	
	in_tval=float(sys.argv[1])
	print res.predict(in_tval), pearsonr(in_arr[:,0], in_arr[:,2])
	
	return 0
	
	#infile, if from command line
	vals=np.loadtxt('../tables/u_flags.txt', usecols=(1, 2), skiprows=1)
	out=[]
	for i in vals:
		arr=ni_ddc(i[0])
		out.append([np.mean(arr), np.std(arr)])
	out=np.array(out)
	print out
	return 0
	#load the nickel mass estimates and the calculated bolometric luminosity
	rr=np.loadtxt("../tables/ni_dif.tex", dtype='string', usecols=(0, 1,2, 3), delimiter='&')
	
	
	
	#load Dm15 for variable rise calculation 
	dm=np.loadtxt(pt+"tmax_dm15.dat", dtype='string')
	
	#change name formats (remove the ' ' at the end)
	nm=[i[:-1] for i in rr[:,0]]; nm=np.array(nm)
	rr[:,0]=nm
	
	#create the arrays for the estimation using the err wala function re
	dmsamp=[[float(i[3]), float(i[4])] for i in dm if i[0] in rr[:,0]]
	lbol=[[float(rr[rr[:,0]==i[0]][0][1]), float(rr[rr[:,0]==i[0]][0][2])] for i in dm if i[0] in rr[:,0]]
	
	
	dmsamp=np.array(dmsamp); lbol=np.array(lbol)
	#coeffcients presuming alpha=1
	factot=[rise_err(dmsamp[i], lbol[i]) for i in range(len(dmsamp))]; factot=np.array(factot)
	
	#calculate the Ni masses
	z=lbol[:,0]/factot[:,0]
	
	#calculate the errors
	zerr=z*((lbol[:,1]/lbol[:,0])+(factot[:,1]/factot[:,0]))
	
	
	print  np.mean(factot[:,1]), np.mean(factot[:,0])
	
	#stop here for the error analysis
	return 0
	
	
	#does bolometric peak calculations 
	
	
main()
