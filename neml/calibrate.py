import numpy as np
import scipy.optimize as opt
import drivers

import scipy.interpolate as inter
import scipy.integrate as qud

import multiprocessing

from lxml import etree as ET

def expand(data, conv = float):
  """
    Expand a string containing either a scalar or a comma separated list
    into either a scalar or a list, hitting each entry with float

    Parameters:
      data      string containing a comma separated list

    Optional:
      conv      converter to use on the string
  """
  spl = data.split(',')
  mspl = map(conv, spl)
  if len(mspl) == 1:
    return mspl[0]
  else:
    return mspl

class Evaluator(object):
  def __init__(self, model_maker, params, weights, penalty):
    self.model_maker = model_maker
    self.params = params
    self.weights = weights
    self.penalty = penalty

  def __call__(self, scase):
    model = self.model_maker(self.params)
    case = ET.fromstring(scase)
    try:
      v = evaluate_case(case, model, weights = self.weights)
    except Exception as e:
      return self.penalty
    if np.isnan(v) or np.isinf(v):
      return self.penalty
    else:
      return v

def evaluate_cases(cases, model_maker, params, weights = {"uniaxial": 1.0,
  "relax_strain": 1.0, "relax_stress": 1.0, "cyclic_stress": 1.0,
  "cyclic_strain": 1.0}, nthreads = 1, penalty = 1.0e4):
  """
    Evaluate a bunch of cases

    Parameters:
      cases     list of xedl cases
      model     model

    Optional:
      weights   weights for each type of experiment
      nthreads  number of threads to use
  """
  # Sigh, lxml objects aren't pickleable.  Like nearly everything...
  string_cases = [ET.tostring(c) for c in cases]

  evaler = Evaluator(model_maker, params, weights, penalty)

  pool = multiprocessing.Pool(nthreads)

  res = list(pool.imap(evaler, string_cases))

  pool.close()
  pool.terminate()
  pool.join()

  return sum(res)

def evaluate_case(case, model, weights = {"uniaxial": 1.0,
  "relax_strain": 1.0, "relax_stress": 1.0, "cyclic_stress": 1.0,
  "cyclic_strain": 1.0}):
  """
    Evaluate a model against an xedl experimental case
    
    Right now supports uniaxial, creep, stress relaxation, and stress/strain
    controlled cyclic tests

    Parameters:
      case      xedl case
      model     model
    
    Optional:
      weights   weights for each type of experiment
  """
  # Determine the type
  if case.attrib['type'] == 'monotonic' and case.find('control').text == 'strain':
    return weights['uniaxial'] * evaluate_uniaxial(case, model)
  elif case.attrib['type'] == 'relaxation' and case.find('control').text == 'strain':
    return weights['relax_strain'] * evaluate_stressrelax(case, model)
  elif case.attrib['type'] == 'relaxation' and case.find('control').text == 'stress':
    return weights['relax_stress'] * evaluate_creep(case, model)
  elif case.attrib['type'] == 'cyclic' and case.find('control').text == 'strain':
    return weights['cyclic_strain'] * evaluate_cyclic_strain(case, model)
  elif case.attrib['type'] == 'cyclic' and case.find('control').text == 'stress':
    return weights['cyclic_stress'] * evaluate_cyclic_stress(case, model)
  else:
    raise ValueError("Unknown experiment with type %s." % case.attrib['type'])

def evaluate_uniaxial(case, model):
  """
    Evaluate a uniaxial test
  """
  T = float(case.find('temperature').text)
  erate = float(case.find('rate').text)
  exp_strain = expand(case.find('strain').text)
  exp_stress = expand(case.find('stress').text)

  emax = np.max(exp_strain)

  res = drivers.uniaxial_test(model, erate, T = T, emax = emax)

  exp_fn = inter.interp1d(exp_strain, exp_stress)
  mod_fn = inter.interp1d(res['strain'], res['stress'])
  
  vs = qud.quad(lambda e: np.abs(exp_fn(e)), 0.0, emax, epsabs = 1.0e-4, epsrel = 1.0e-2,
      full_output = 1)
  v1 = vs[0]
  vs = qud.quad(lambda e: np.abs(exp_fn(e) - mod_fn(e)), 
      0.0, emax, epsabs = 1.0e-4, epsrel = 1.0e-2, full_output = 1)
  v2 = vs[0]

  return v2 / v1

def cost_uniaxial(case):
  return 250

def evaluate_stressrelax(case, model):
  """
    Evaluate a stress relaxation test
  """
  T = float(case.find('temperature').text)
  erate = float(case.find('rate').text)
  strain = float(case.find('value').text)
  
  exp_time = expand(case.find('time').text)
  exp_relax = expand(case.find('relax').text)
  
  min_time = np.min(exp_time)
  max_time = np.max(exp_time)

  res = drivers.stress_relaxation(model, strain, erate, max_time, 
      T = T)
  
  exp_fn = inter.interp1d(exp_time, exp_relax)
  mod_fn = inter.interp1d(res['rtime'], res['rstress'])

  vs = qud.quad(lambda e: np.abs(exp_fn(e)), min_time, max_time*0.99, epsabs = 1.0e-4, 
      epsrel = 1.0e-3, full_output = 1)
  v1 = vs[0]
  vs = qud.quad(lambda e: np.abs(exp_fn(e) - mod_fn(e)), 
      min_time*1.01, max_time*0.99, epsabs = 1.0e-4, epsrel = 1.0e-3,
      full_output = 1)
  v2 = vs[0]

  return v2 / v1

def cost_stressrelax(model):
  return 900

def evaluate_creep(case, model):
  """
    Evaluate a creep test
  """
  T = float(case.find('temperature').text)
  srate = float(case.find('rate').text)
  stress = float(case.find('value').text)

  exp_time = expand(case.find('time').text)
  exp_strain = expand(case.find('relax').text)
  
  min_time = np.min(exp_time)
  max_time = np.max(exp_time)

  res = drivers.creep(model, stress, srate, max_time, T = T)

  if len(res['rtime']) < 2:
    raise Exception('No good')

  max_time = min([max_time, np.max(res['rtime'])])
  min_time = max([min_time, np.min(res['rtime'])])

  exp_fn = inter.interp1d(exp_time, exp_strain)
  mod_fn = inter.interp1d(res['rtime'], res['rstrain'])
  
  vs = qud.quad(lambda e: np.abs(exp_fn(e)), min_time*1.01, max_time*0.99, epsabs = 1.0e-4, 
      epsrel = 1.0e-3, full_output = 1)
  v1 = vs[0]
  vs = qud.quad(lambda e: np.abs(exp_fn(e) - mod_fn(e)), 
      min_time*1.01, max_time*0.99, epsabs = 1.0e-4, epsrel = 1.0e-3,
      full_output = 1)
  v2 = vs[0]

  return v2/v1

def cost_creep(case):
  return 900

def evaluate_cyclic_strain(case, model):
  """
    Evaluate a strain controlled cyclic test
  """
  T = float(case.find('temperature').text)
  rate = float(case.find('rate').text)
  strain = float(case.find('value').text)
  ratio = float(case.find('ratio').text)

  exp_cycle = expand(case.find('cycle').text)
  if case.find('max') is None:
    return 0.0
    #raise RuntimeError("Implement this please...")
  exp_max = expand(case.find('max').text)

  if case.find('hold') is not None:
    #raise RuntimeError("Holds on strain controlled tests not implemented...")
    return 0
  
  max_cycles = np.max(exp_cycle)
  min_cycles = np.min(exp_cycle)

  # This is a hack, please fix for the love of god
  if max_cycles > 2000:
    max_cycles = 2000

  res = drivers.strain_cyclic(model, strain, ratio, rate,
      int(max_cycles)+1, T = T)

  max_cycles = min([max_cycles, np.max(res['cycles'])])
  min_cycles = max([min_cycles, np.min(res['cycles'])])

  exp_fn = inter.interp1d(exp_cycle, exp_max)
  mod_fn = inter.interp1d(res['cycles'], res['max'])

  vs = qud.quad(lambda e: np.abs(exp_fn(e)), min_cycles, int(max_cycles)*0.99, epsabs = 1.0e-4, 
      epsrel = 1.0e-3, full_output = 1)
  v1 = vs[0]
  vs = qud.quad(lambda e: np.abs(exp_fn(e) - mod_fn(e)), 
      min_cycles, int(max_cycles)*0.99, epsabs = 1.0e-4, epsrel = 1.0e-3,
      full_output = 1)
  v2 = vs[0]

  return v2/v1

def cost_cycle_strain(case):
  nsteps = np.max(expand(case.find('cycle').text))
  return 50 * nsteps

def evaluate_cyclic_stress(case, model):
  """
    Evaluate a stress controlled cyclic test
  """
  T = float(case.find('temperature').text)
  rate = float(case.find('rate').text)
  stress = float(case.find('value').text)
  ratio = float(case.find('ratio').text)

  exp_cycle = expand(case.find('cycle').text)

  if case.find('max') is not None:
    compare = 'max'
  elif case.find('mean') is not None:
    compare = 'mean'
  elif case.find('min') is not None:
    compare = 'min'
  else:
    return 0.0
  exp_data = expand(case.find(compare).text)

  if case.find('hold') is not None:
    if case.find('hold').attrib['type'] == "max":
      hold = [float(case.find('hold').text), 0.0]
    elif case.find('hold').attrib['type'] == "min":
      hold = [0.0, float(case.find('hold').text)]
    else:
      raise ValueError("Unknown hold type!")
  else:
    hold = None
  
  max_cycles = np.max(exp_cycle)
  min_cycles = np.min(exp_cycle)

  # This is a hack, please fix for the love of god
  if max_cycles > 2000:
    max_cycles = 2000

  res = drivers.stress_cyclic(model, stress, ratio, rate, int(max_cycles)+1,
      T = T, hold_time = hold)

  exp_fn = inter.interp1d(exp_cycle, exp_data)
  mod_fn = inter.interp1d(res['cycles'], res[compare])
  
  vs = qud.quad(lambda e: np.abs(exp_fn(e)), min_cycles, max_cycles*0.99, epsabs = 1.0e-4, 
      epsrel = 1.0e-3, full_output = 1)
  v1 = vs[0]
  vs = qud.quad(lambda e: np.abs(exp_fn(e) - mod_fn(e)), 
      min_cycles, max_cycles*0.99, epsabs = 1.0e-4, epsrel = 1.0e-3,
      full_output = 1)
  v2 = vs[0]

  return v2/v1

def cost_cycle_stress(case):
  nsteps = np.max(expand(case.find('cycle').text))
  return 50 * nsteps
