#!/usr/bin/env python

import numpy as np

import drivers

def yanli_cooling_cycle(T1, T2, dT, hold, delay):
  """
    Generate Yanli's cycle with the delay on the cooling side

    Parameters:
      T1        cool bar temperatures
      T2        hot bar temperatures
      dT        heating/cooling rate
      hold      hold time at same temperature
      delay     delay between cooling
  """
  httime = np.abs(T2-T1) / np.abs(dT)
  per = 2.0 * httime + hold + delay
  
  T1f = lambda t: np.piecewise(t,
      [t<=httime, np.logical_and(t<=(httime+hold), t>httime), 
        np.logical_and(t<=(2.0*httime+hold), t>(httime+hold)),
          t>(2.0*httime+hold)],
      [lambda tt: tt*dT + T1, T2, lambda tt: T2 - dT*(t-httime-hold), T1])

  T2f = lambda t: np.piecewise(t,
      [t<=httime, np.logical_and(t<=(httime+hold), t>httime), 
        np.logical_and(t<=(httime+hold+delay), t>(httime+hold)),
          t>(httime+hold+delay)],
      [lambda tt: tt*dT + T1, T2, T2, lambda tt: T2 - dT*(t-httime-hold-delay)])

  return T1f, T2f, per

def yanli_heating_cycle(T1, T2, dT, hold, delay):
  """
    Generate Yanli's cycle with the delay on the heating side

    Parameters:
      T1        cool bar temperatures
      T2        hot bar temperatures
      dT        heating/cooling rate
      hold      hold time at same temperature
      delay     delay between heating
  """
  httime = np.abs(T2-T1) / np.abs(dT)
  per = 2.0 * httime + hold + delay
  
  T1f = lambda t: np.piecewise(t,
      [t<=httime, np.logical_and(t>httime, t<=(httime+delay+hold)), 
        t>(httime+delay+hold)],
      [lambda tt: tt*dT + T1, T2, lambda tt: T2 - dT*(tt-httime-hold-delay)])

  T2f = lambda t: np.piecewise(t,
      [t<=delay, np.logical_and(t<=(delay+httime), t>delay), 
        np.logical_and(t<=(delay+httime+hold), t>(delay+httime)),
          t>(httime+hold+delay)],
      [T1, lambda tt: (tt-delay)*dT + T1, T2, lambda tt: T2 - dT*(t-httime-hold-delay)])

  return T1f, T2f, per

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

    return True, res

  else:
    return False, res

