#!/usr/bin/env python

import numpy as np

import drivers

def epp_two_bar(models, A1, A2, T1, T2, period, P, load_time, ncycles):
  """
    Run the EPP procedure over a two bar system

    Parameters:
      models    dictionary of x: model
      A1        area of bar 1
      A2        area of bar 2
      T1        temperature cycle bar 1
      T2        temperature cycle bar 2
      period    period of load cycle
      P         common load
      load_time ramp time for load
      ncycles   number of cycles to run
  """
  for x, model in models.iteritems():
    res = drivers.twobar_test(model, A1, A2, T1, T2, period, P, load_time,
        ncycles)
    if res['classification'] not in ('elastic', 'elastic shakedown', 
        'plastic shakedown'):
      continue

    strains = (res['strain_inelastic_1'], res['strain_inelastic_2'])
    big = np.max(strains)
    sma = np.min(strains)

    if (big + x > 0.05) or (sma + x > 0.01):
      continue

    return True

  else:
    return False

