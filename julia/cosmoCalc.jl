H0=70
WM=0.27
WV=0.73
WK=1-WM-WV
WR=0
c = 299792.458

function mod(z)
	az = 1.0/(1+1.0*z)
	age = 0
	n=1000        
	

	DTT = 0.0
	DCMR = 0.0
	for i = 1:n
		a = az+(1-az)*(i+0.5)/n
  		adot = sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))
 		DTT = DTT + 1./adot
  		DCMR = DCMR + 1./(a*adot)	
	end
	DTT = (1.-az)*DTT/n
	DCMR = (1.-az)*DCMR/n

	DCMR_Mpc = (c/H0)*DCMR
	ratio = 1.00
	x = sqrt(abs(WK))*DCMR

	if x >0.1
		if WK > 0
			ratio =  0.5*(exp(x)-exp(-x))/x
		else

			ratio= sin(x)/x		
		end	

	else
		y=x*x
	end	
	if  WK < 0
		y=-y
	end
	ratio = 1. + y/6. + y*y/120.
	DCMT = ratio*DCMR
	DA = az*DCMT	
	DL = DA/(az*az)
	DL_Mpc = (c/H0)*DL
	mu=5*log10(DL_Mpc)+25
	return mu
	
end
