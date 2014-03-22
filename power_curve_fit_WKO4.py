import pylab as pl
import numpy as np
from scipy import optimize
from scipy.stats.distributions import  t

#http://kitchingroup.cheme.cmu.edu/blog/2013/02/12/Nonlinear-curve-fitting-with-parameter-confidence-intervals/
#http://nbviewer.ipython.org/url/media.usm.maine.edu/~pauln/ScipyScriptRepo/CurveFitting.ipynb

ttimes= [  '1s', '8s', '15s', '1m', '1m30','2m','3m','4m','5m','6m','10m','20m','30m','45m','1h','1h30','2h', '3h', '4h', '5h', '6h',' 7h', '8h']
times = [   1.0,     8, 15.0, 60.0,     90, 120, 180, 240, 300, 360,  600, 1200, 1800, 2700,3600,  5400,7200,10800,14400,18000,21600,25200,28800]
power = [1123.0,   842,775.0,639.0,    535, 409, 401, 371, 343, 327,  320,  299,  287,  277, 270,   239, 224,  180,  173,  168,  167,  165,  159]

pl.xscale('log')
pl.scatter(times,power,c='k',label='data')
pl.xlabel('Duration')
pl.ylabel('Power (Watt)')
pl.xticks(times, ttimes)
pl.xlim([1, 30000])
pl.ylim([0,1200])
pl.grid(True)
pl.title('Power duration', fontsize=14, fontweight='bold')


def fpower(t,Pow1, tau1,Pow2,tau2):
	return Pow1*tau1/(t+tau1) + Pow2*tau2/(t+tau2)

# WKO4 >?? : http://veloclinic.tumblr.com/post/78082056106/hey-cog-is-the-wko4-model-this
def fpower1(t, FRC, Pmax, FTP, alpha):
	return FRC/t*(1-np.exp(-t/(FRC/(Pmax-FTP)))) + FTP + alpha*(t-3600)

X = np.array(times)
Y = np.array(power)
popt, pcov= optimize.curve_fit(fpower1, X, Y, [300.0, 1123.0, 280.0, 0.9])
FRC, Pmax, FTP, alpha = popt
Y_fit = map(lambda t: fpower1(t,FRC, Pmax, FTP, alpha) , times)
pl.plot(times, Y_fit, c='g', label='fit', linewidth=2.5)


alphac = 0.05 # 95% confidence interval = 100*(1-alpha)
n = len(power)    # number of data points
p = len(popt) # number of parameters
dof = max(0, n - p) # number of degrees of freedom
# student-t value for the dof and confidence level
tval = t.ppf(1.0-alphac/2., dof)

for i, p,var in zip(range(n), popt, np.diag(pcov)):
    sigma = var**0.5
    print 'p{0}: {1} [{2}  {3}]'.format(i, p,
                                  p - sigma*tval,
                                  p + sigma*tval)
print 'FRC =', FRC, ' Watt'
print 'Pmax =', Pmax, ' sec'
print 'FTP =', FTP, ' Watt'
print 'alpha =', alpha, ' sec'

lower = map(lambda p,var: p - (var**0.5)*tval , popt, np.diag(pcov))
upper = map(lambda p,var: p + (var**0.5)*tval , popt, np.diag(pcov))
yfit = map(lambda t: fpower1(t, *lower) , times)
pl.plot(times,yfit,'--', color='c')
yfit = map(lambda t: fpower1(t, *upper) , times)
pl.plot(times,yfit,'--', color='c', label='CI 95%')
pl.legend(['fit','CI 95%'],loc='upper right')


figure = pl.gcf() # get current figure
figure.set_size_inches(18, 6)
# when saving, specify the DPI
pl.savefig("power_WKO4.png", dpi = 200)
#pl.savefig("power.svg", dpi = 200)
pl.show()






