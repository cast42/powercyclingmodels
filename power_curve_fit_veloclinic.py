import pylab as pl
import numpy as np
from scipy import optimize
from scipy.stats.distributions import  t

#http://kitchingroup.cheme.cmu.edu/blog/2013/02/12/Nonlinear-curve-fitting-with-parameter-confidence-intervals/

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

# model from veloclinic: https://docs.google.com/document/d/1-JJWsO-xJhTjhX5mn-IweTMIb0i4zjCwuQAPsnK0OBc/edit
def fpower(t,Pow1, tau1,Pow2,tau2):
	return Pow1*tau1/(t+tau1) + Pow2*tau2/(t+tau2)

X = np.array(times)
Y = np.array(power)
popt, pcov= optimize.curve_fit(fpower, X, Y)
Pow2, tau2, Pow1, tau1 = popt
Y_fit = map(lambda t: fpower(t,Pow1,tau1,Pow2,tau2) , times)
pl.plot(times, Y_fit, c='g', label='fit', linewidth=2.5)


alpha = 0.05 # 95% confidence interval = 100*(1-alpha)
n = len(power)    # number of data points
p = len(popt) # number of parameters
dof = max(0, n - p) # number of degrees of freedom
# student-t value for the dof and confidence level
tval = t.ppf(1.0-alpha/2., dof)

for i, p,var in zip(range(n), popt, np.diag(pcov)):
    sigma = var**0.5
    print 'p{0}: {1} [{2}  {3}]'.format(i, p,
                                  p - sigma*tval,
                                  p + sigma*tval)
print 'Pow1 =', Pow1, ' Watt'
print 'tau1 =', tau1, ' sec'
print 'Pow2 =', Pow2, ' Watt'
print 'tau2 =', tau2, ' sec'

lower = map(lambda p,var: p - (var**0.5)*tval , popt, np.diag(pcov))
upper = map(lambda p,var: p + (var**0.5)*tval , popt, np.diag(pcov))
yfit = map(lambda t: fpower(t, *lower) , times)
pl.plot(times,yfit,'--', color='c')
yfit = map(lambda t: fpower(t, *upper) , times)
pl.plot(times,yfit,'--', color='c', label='CI 95%')
pl.legend(['fit','CI 95%'],loc='upper right')


figure = pl.gcf() # get current figure
figure.set_size_inches(18, 6)
# when saving, specify the DPI
pl.savefig("power_veloclinic.png", dpi = 200)
#pl.savefig("power.svg", dpi = 200)
pl.show()






