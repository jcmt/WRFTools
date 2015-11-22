

def scramble(orig):
  '''
  Same as random.shuffle from Python but it keeps the original variable
  '''
  import random
  import numpy as np

  dest = np.copy(orig)
  random.shuffle(dest)
  return dest

def erroreval(obs, model):
  import numpy as np
  import pandas as pd

  ErrStat = {}

  thetal = model - obs
  thetal2 = thetal ** 2

  theta0 = model[~np.isnan(model)].mean()
  tmodel0 = model - theta0
  theta0_obs = obs[~np.isnan(obs)].mean()
  tobs0 = obs - theta0_obs
  difft0 = (tmodel0 - tobs0) ** 2

  ErrStat['BIAS'] = thetal[~np.isnan(thetal)].mean()
  ErrStat['E'] = thetal2[~np.isnan(thetal2)].mean() ** 0.5
  ErrStat['STDE'] = thetal[~np.isnan(thetal)].std()
  ErrStat['E_UB'] = difft0[~np.isnan(difft0)].mean() ** 0.5
  ErrStat['S'] = tmodel0[~np.isnan(tmodel0)].std()
  ErrStat['S_obs'] = tobs0[~np.isnan(tobs0)].std()

  model = pd.Series(tmodel0)
  obs = pd.Series(tobs0)
  ErrStat['R'] = model.corr(obs)

  return ErrStat


################# MangaCbar improved by MMA #######################
def read_file(f):
  import numpy as np

  l = open(f).readlines()

  i = -1
  for k in l:
    i += 1
    if k.startswith('static'):
      break

  l = l[i:]

  i0 = l[0].find('{')
  i1 = l[-1].find('}')

  l[0] = l[0][i0+1:]
  l[-1] = l[-1][:i1]

  r = []
  for i in range(len(l)):
    line = l[i].replace(',', ' ').strip()
    vals = [int(j) for j in line.split()]
    r = np.hstack((r, vals))

  r.shape = r.size/3, 3
  return r/255.


def gencmap(file, name='cmap', N=None):
  import pylab as pl
  import os

  path = os.path.join(os.path.join(os.path.dirname(__file__)), 'colormaps/')
  r = read_file(path + file)
  return pl.matplotlib.colors.ListedColormap(r, name=name, N=N)


def smooth(x, window_len=11, window='hanning'):
  """
  smooth the data using a window with requested size.

  This method is based on the convolution of a scaled window with the signal.
  The signal is prepared by introducing reflected copies of the signal
  (with the window size) in both ends so that transient parts are minimized
  in the begining and end part of the output signal.

  input:
      x: the input signal
      window_len: the dimension of the smoothing window; should be an odd integer
      window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
          flat window will produce a moving average smoothing.

  output:
      the smoothed signal

  example:

  t=linspace(-2,2,0.1)
  x=sin(t)+randn(len(t))*0.1
  y=smooth(x)

  see also:

  numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
  scipy.signal.lfilter

  TODO: the window parameter could be the window itself if an array instead of a string
  NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
  """
  import numpy as np

  if x.ndim != 1:
    raise(ValueError, "smooth only accepts 1 dimension arrays.")

  if x.size < window_len:
    raise(ValueError, "Input vector needs to be bigger than window size.")

  if window_len < 3:
    return x

  if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
    raise(ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

  s = np.r_[x[window_len-1:0:-1], x, x[-1:-window_len:-1]]
  #print(len(s))
  if window == 'flat': #moving average
    w = np.ones(window_len, 'd')
  else:
    w = eval('np.'+window+'(window_len)')

  y = np.convolve(w/w.sum(), s, mode='valid')
  CutTail = int(np.fix(window_len/2.0))
  y = y[CutTail:-CutTail]

  return y
