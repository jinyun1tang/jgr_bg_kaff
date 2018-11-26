#!/usr/bin/env python

"""python script for linearfitting between x and y
Author: Jinyun Tang <jinyuntang@gmail.com>
"""

import numpy as np
import math
from numpy.matlib import rand,zeros,ones,empty,eye

def linearfit(x,y):
    fit = np.polyfit(x,y,1)
    fit_fn = np.poly1d(fit)
    w1=np.zeros(6)
    w=np.zeros(8)
    w1[0]=fit[0]
    w1[1]=fit[1]

    yhat=fit_fn(x)
    ybar=np.mean(yhat)
    sst=np.dot(y-ybar,y-ybar)
    sse=np.dot(ybar-yhat,ybar-yhat)
    ssr=np.dot(y-yhat,y-yhat)
    n=len(x)
    w1[2]=np.sqrt(np.dot(x-y,x-y)/n)
    w1[3]=sse/sst
    "do F-test"
    mse=sse
    msr=ssr/(n-2.0)
    ftest=mse/msr
    from scipy.stats import f
    p = f.cdf(ftest, 1, n - 2);
    w1[4]=1.0-p
    sgm1=np.sqrt((np.dot(y,y)-w1[1]*np.sum(y)-w1[0]*np.dot(x,y))/(n-2.0))
    sgm2=np.sqrt(n*np.dot(x,x)-np.sum(x)*np.sum(x))
    sgm=sgm1*sgm2/(sgm2*sgm2+1.e-20)
    "intercept"
    sgb=sgm*np.sqrt(np.dot(x,x)/n)
    xa=x-np.mean(x)
    ya=y-np.mean(y)
    w1[5]=np.sum(np.dot(xa,ya))/np.sqrt(np.dot(xa,xa))/np.sqrt(np.dot(ya,ya))
    w[0]=w1[0]  #slope
    w[1]=sgm    #standard dev of slope
    w[2]=w1[1]  #intercept
    w[3]=sgb    #standard dev of intercept
    w[4]=w1[2]  #root mean square error
    w[5]=w1[3]  #R-squared
    w[6]=w1[4]  #p value
    w[7]=w1[5]  #correlation coefficient
    return w

""" example usage
x = np.asarray([1,2,3,4])
y = np.asarray([3,5,7,9]) # 10, not 9, so the fit isn't perfect

fit=linearfit(x,y)

print(fit)
"""
