'''
Find the centroid 
Written by MJ
'''
from __future__ import division
import numpy as np
from scipy import signal , linalg
from scipy.linalg import cho_factor, cho_solve
from sep import extract

x, y = np.meshgrid(range(-1, 2), range(-1, 2), indexing="ij")
x, y = x.flatten(), y.flatten()
AT = np.vstack((x*x, y*y, x*y, x, y, np.ones_like(x)))

C = np.identity(9)
ATA = np.dot(AT, np.dot(np.linalg.inv(C) , AT.T))
factor = cho_factor(ATA, overwrite_a=True)
#ATA = np.dot(AT, AT.T)
#factor = cho_factor(ATA, overwrite_a=True)


def fit_3x3(im):
    imgg = np.dot(AT , np.dot(np.linalg.inv(C) , im.flatten()))
    a, b, c, d, e, f = cho_solve(factor, imgg)
    m = 1. / (4 * a * b - c*c)
    x = (c * e - 2 * b * d) * m
    y = (c * d - 2 * a * e) * m
    return x, y

def makeGaussian(size , FWHM , e = 0 , center = None):

    f = FWHM/(2.35482)
    x = np.linspace(0.5, size-.5 , size)
    y = x[:,np.newaxis]

    if center is None:
        x0 = y0 = size / 2.
    else:
        x0 = center[0]
        y0 = center[1]
        
    r = ((x-x0)**2. + ((y-y0)**2.)*(1. + np.abs(e))**2./(1. - np.abs(e))**2.)**.5
    factor = 1./(2.*np.pi*f**2.)   
    return factor*np.exp((-1.*r**2.)/(2*f**2.))


def MAD(a, axis=None):
    """Compute the median absolute deviation"""
    a = np.array(a, copy=False)
    a_median = np.median(a, axis=axis)

    #re-broadcast the output median array to subtract it
    if axis is not None:
        shape = list(a_median.shape)
        shape.append(1)
        a_median = a_median.reshape(shape)

    #calculated the median average deviation
    return np.median(np.abs(a - a_median), axis=axis)/0.6745

def find_centroid(data, t):
  '''
  filter_kernel = makeGaussian(17, 5. , 0 , np.array([8.5,8.5]))
  source = extract(data , t, filter_kernel=filter_kernel)
  '''
  #t=20
  #source = extract(data , t)
  source = extract(data, 100)
  #source = extract(data, t)
  #source = extract(data, 1)
 
  a = source['a']
  b = source['b']
  print 'a: {0}'.format(a)
  print 'b: {0}'.format(b)
  flux = source['cflux']
  #arg = np.argsort(flux)[-1]
  try:
    arg = np.argsort(flux)[-1]
  except IndexError:
    return None
  print flux

  try:
    #fwhm = np.sqrt(np.max(a)*np.max(b))
    fwhm = np.sqrt(a[arg]*b[arg])
    print 'fwhm:{0}'.format(fwhm)
  except ValueError:
    return None
  size = data.shape[0]
  zero = size/2 + .5
  kernel = makeGaussian(17, fwhm , 0 , np.array([8.5,8.5]))
  img = signal.convolve2d(data , kernel , mode = "same")
  max_value = np.max(img)
  xi, yi = np.unravel_index(np.argmax(img), img.shape)
  
  if (xi >= 1 and xi < img.shape[0] - 1 and yi >= 1 and yi < img.shape[1] - 1):
      ox, oy = fit_3x3(img[xi-1:xi+2, yi-1:yi+2])
  else:
      ox , oy = 0. , 0.
  #return xi + ox + .5 , yi + oy + .5, max_value, flux[arg[-1]]
  return xi + ox + .5 , yi + oy + .5, max_value, flux[arg]

def find_source(data, t):
  source = extract(data , t)
  a = source['a']
  b = source['b']
  x = source['x']
  y = source['y']
  flux = source['flux']
  print x
  print y
  print flux
  arg = np.argsort(flux)
  print x[arg[-1]], y[arg[-1]], flux[arg[-1]]
  try:
    fwhm = np.sqrt(np.max(a)*np.max(b))
    print fwhm
  except ValueError:
    return None

  return x[arg[-1]], y[arg[-1]], flux[arg[-1]]


if __name__ == "__main__":
    print 'c3 main'
