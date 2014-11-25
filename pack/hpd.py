## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.
#

""" Calculate Bayesian statistics sucj as Heighst Posterior Density (HPD) and
Effective Sample Size (Beast interpretation).
"""

from __future__ import division

import numpy

__all__ = ["hpd", "effectiveSampleSize"]

def hpd(data, level=0.68) :
  """ The Highest Posterior Density (credible) interval of C{data} at level C{level}
  (0 < level < 1). """ 
  
  d = list(data)
  d.sort()

  nData = len(data)
  nIn = int(round(level * nData))
  if nIn < 2 :
    raise RuntimeError("not enough data")
  
  i = 0
  r = d[i+nIn-1] - d[i]
  for k in range(len(d) - (nIn - 1)) :
    rk = d[k+nIn-1] - d[k]
    if rk < r :
      r = rk
      i = k

  assert 0 <= i <= i+nIn-1 < len(d)
  
  return (d[i], d[i+nIn-1])
