module dist_cosm

export lum_dist
c=2.99e5
function lum_dist( z, om, ol, h)
	q0=(om/2)-ol
	a=z*(1-q0)/(np.sqrt(1+2*q0*z)+1+q0*z)
	dl=(1+a)*c*z/h0
	return dl/(1+z)
end


end #module
