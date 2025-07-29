import numpy as np
from scipy.stats import genlogistic, gompertz
import matplotlib.pyplot as plt

from get import insult_severity



def get_hist(InsultDistribution):
       
       CONFIG = dict(
              InsultDistribution = InsultDistribution,
              LogorithmicPopulationCap = False
       )
       fig, ax = plt.subplots(1, 1)

       r = [insult_severity(CONFIG) for i in range(1000)]

       plt.title(InsultDistribution)
       ax.set_xlim([0,150])

       ax.hist(r, density=True, bins='auto', histtype='barstacked', log=True, cumulative=False)
       plt.savefig('graphs/150'+InsultDistribution+ 'Histogram.png')


# get_hist('T3')     
# get_hist('T4')     
# get_hist('T5')     
# get_hist('B2')     

  

def plot (c, loc1, cdf_plt=False, pdf_plt=False, use_abs=False, use_log_scale=False, fcn=genlogistic):

       
       # x = np.linspace(0,
       #                 10, 100) 
       # x = np.linspace(fcn.ppf(0.01, c, loc = loc1),
       #               fcn.ppf(0.99, c, loc = loc1), 1000) 
       
       x = np.linspace(0,
                     fcn.ppf(0.99, c, loc = loc1), 100) 
       if use_abs:
              b = [fcn.pdf(-i, c, loc = loc1) + fcn.pdf(i, c, loc = loc1)  for i in x]

       if pdf_plt:
              fig, ax = plt.subplots(1, 1)

              mean, var, skew, kurt = fcn.stats(c, moments='mvsk', loc = loc1)
              # print(mean, var, skew, kurt )
              # Display the cumulative density function (cdf):

              if use_abs:
                     ax.plot(x, b,
                            'r-', lw=5, alpha=0.6, label='fcn pdf')
                            # Alternatively, the distribution object can be called (as a function) to fix the shape, location and scale parameters. This returns a “frozen” RV object holding the given parameters fixed.

                     # Freeze the distribution and display the frozen pdf:

                     rv = fcn(c, loc = loc1)
                     ax.plot(x, b, 'k-', lw=2, label='frozen pdf')
              else :
                     ax.plot(x, fcn.pdf(x, c, loc = loc1),
                            'r-', lw=5, alpha=0.6, label='fcn pdf')
                     # Alternatively, the distribution object can be called (as a function) to fix the shape, location and scale parameters. This returns a “frozen” RV object holding the given parameters fixed.

                     # Freeze the distribution and display the frozen pdf:

                     rv = fcn(c, loc = loc1)
                     ax.plot(x, rv.pdf(x), 'k-', lw=2, label='frozen pdf')

              r = fcn.rvs(c, size=1000, loc = loc1)
              if use_abs:
                     r = [abs(i) for i in r]

              # And compare the histogram:

              ax.hist(r, density=True, bins='auto', cumulative=False, histtype='stepfilled', alpha=0.2)
              ax.set_xlim([x[0], x[-1]])

              # ax.legend(loc='best', frameon=True)
              if use_log_scale:
                     plt.yscale('log')
                     
              t = 'PDF of c: ' + str(c) +  ', loc: ' + str(loc1) + 'with curve ' + str(fcn.name)
              t += ' Using Abs' if use_abs else ''
              plt.title(t)
       
       if cdf_plt:
              fig, ax = plt.subplots(1, 1)

              mean, var, skew, kurt = fcn.stats(c, moments='mvsk')
              # print(mean, var, skew, kurt )
              # fcn.
              # Display the cumulative density function (cdf):
              # loc1 = 10

              if use_abs:
                     b2 = [sum(b[:i])/x[-1] for i in range(len(b))]
                     ax.plot(x, b2,
                            'r-', lw=5, alpha=0.6, label='fcn cdf')
                                   # Alternatively, the distribution object can be called (as a function) to fix the shape, location and scale parameters. This returns a “frozen” RV object holding the given parameters fixed.

                     # Freeze the distribution and display the frozen cdf:

                     rv = fcn(c, loc = loc1)
                     ax.plot(x, b2, 'k-', lw=2, label='frozen cdf')
              else :
                     ax.plot(x, fcn.cdf(x, c, loc = loc1),
                            'r-', lw=5, alpha=0.6, label='fcn cdf')
                     # Alternatively, the distribution object can be called (as a function) to fix the shape, location and scale parameters. This returns a “frozen” RV object holding the given parameters fixed.

                     # Freeze the distribution and display the frozen cdf:

                     rv = fcn(c, loc = loc1)
                     ax.plot(x, rv.cdf(x), 'k-', lw=2, label='frozen cdf')
              # Check accuracy of cdf and ppf:

              # Generate random numbers:

              r = fcn.rvs(c, size=1000, loc = loc1)
              if use_abs:
                     r = [abs(i) for i in r]

              # And compare the histogram:

              ax.hist(r, density=True, bins='auto', cumulative=True, histtype='stepfilled', alpha=0.2)
              # ax.set_xlim([-1, x[-1]])
              ax.set_xlim([x[0], x[-1]])

              # ax.legend(loc='best', frameon=True)

              if use_log_scale:
                     plt.yscale('log')

              t = 'CDF of c: ' + str(c) +  ', loc: ' + str(loc1) + 'with curve ' + str(fcn.name)
              t += ' Using Abs' if use_abs else ''
              plt.title(t)
              #plt.show()


              # Calculate the first four moments:


# plot(.4, 10, cdf_plt=True)
# plot(0.5799,10, cdf_plt=True)
# plot(.7,7, cdf_plt=True)
# plot(.7,8, cdf_plt=True)
# plot(.2,10, pdf_plt=True, cdf_plt=True, use_log_scale=False)
# plot(.5799, 0, pdf_plt=True, cdf_plt=False, use_abs = True)


# for i in range(10,20):
#        plot(i/10, 0, pdf_plt=True, fcn=gompertz)
plt.show()
