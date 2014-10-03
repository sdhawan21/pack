from mag2fl import conv			# regular imports, important one being gp from scikit. mag2fl is a convenient way to read the light curve
from sklearn.gaussian_process import GaussianProcess as gp

import numpy as np
import matplotlib.pyplot as plt
import sys
ng=float(sys.argv[1])		#define nugget and theta0
th=float(sys.argv[2])
band=sys.argv[3]		#input SN and band
obj=sys.argv[4]
lc1=conv().rd_lc(obj, band)
x=lc1['MJD'][:, None]#[lc1['MJD']<54440]
y=lc1[band][:, None]#[lc1['MJD']<54440]
#x=x[1:-3:2]
#y=y[1:-3:2]
g=gp(theta0=th, nugget=ng)
g.fit(x, y)
x_pred=np.linspace(x.min(), x.max())[:,None]
y_pred, MSE=g.predict(x_pred, eval_MSE=True)
print min(y_pred)
plt.plot(x, y, 'b.', label='data')
plt.plot(x_pred, y_pred, 'k:', label='gp fit')
plt.ylim(20, 13)
plt.xlabel('MJD')
plt.ylabel('$m$'+band)
plt.legend(loc=0)
plt.show()
