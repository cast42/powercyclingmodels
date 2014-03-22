import pylab as pl
import numpy as np
from scipy import optimize
from scipy.stats.distributions import  t

#http://kitchingroup.cheme.cmu.edu/blog/2013/02/12/Nonlinear-curve-fitting-with-parameter-confidence-intervals/
#http://nbviewer.ipython.org/url/media.usm.maine.edu/~pauln/ScipyScriptRepo/CurveFitting.ipynb

ttimes= [  '1s', '4s','8s', '15s', '1m', '1m30','2m','3m','4m','5m','6m','10m','20m','30m','45m','1h','1h30','2h', '3h', '4h', '5h', '6h',' 7h', '8h']
times = [   1.0,  4.0,   8, 15.0, 60.0,     90, 120, 180, 240, 300, 360,  600, 1200, 1800, 2700,3600,  5400,7200,10800,14400,18000,21600,25200,28800]
power = [1123.0,961.0, 842,775.0,639.0,    535, 409, 401, 371, 343, 327,  320,  299,  287,  277, 270,   239, 224,  180,  173,  168,  167,  165,  159]

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

# model from WKO4 ???
def fpower1(t, FRC, Pmax, FTP, alpha):
	return FRC/t*(1-np.exp(-t/(FRC/(Pmax-FTP)))) + FTP + alpha*(t-3600)

# model from Djconnel : http://djconnel.blogspot.be/2014/02/power-model-summary.html
def fpower2(t, P1,tau1, P2,tau2, alpha2):
	# P|t = P1 ( tau1 / t) ( 1 - exp[ -t / tau1 ] ) + P2 / ( 1 + t / tau2 ) * alpha2
	return P1*(tau1/t)*(1 - np.exp(-t/tau1)) + P2/np.power((1+ t/tau2),alpha2)

X = np.array(times)
Y = np.array(power)
popt, pcov= optimize.curve_fit(fpower, X, Y)
Pow2, tau2, Pow1, tau1 = popt
Y_fit = map(lambda t: fpower(t,Pow1,tau1,Pow2,tau2) , times)
pl.plot(times, Y_fit, c='g', label='veloclinic', linewidth=2.5)

popt1, pcov1 = optimize.curve_fit(fpower1, X, Y, [300.0, 1123.0, 280.0, 0.9])
FRC, Pmax, FTP, alpha = popt1
Y_fit1 = map(lambda t: fpower1(t,FRC, Pmax, FTP, alpha) , times)
pl.plot(times, Y_fit1, c='r', label='WKO4', linewidth=2.5)

popt2, pcov2 = optimize.curve_fit(fpower2, X, Y, [1200.0, 300.0, 280.0, 25000.0, 0.5])
P1, tau1, P2, tau2, alpha2 = popt2
Y_fit2 = map(lambda t: fpower2(t,P1,tau1, P2,tau2, alpha2) , times)
pl.plot(times, Y_fit2, c='b', label='djconnel', linewidth=2.5)

pl.legend(loc='upper right')

figure = pl.gcf() # get current figure
figure.set_size_inches(18, 6)
# when saving, specify the DPI
pl.savefig("power_all.png", dpi = 200)
#pl.savefig("power.svg", dpi = 200)
pl.show()






