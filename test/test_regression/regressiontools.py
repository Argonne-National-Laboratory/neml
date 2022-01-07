#!/usr/bin/env python3

import os.path
import sys
import argparse

import tempfile

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from neml import models, parse, drivers

default_xml_file = "model.xml"
default_model_name = "model"
default_result_name = "data.csv"

def generate_tensile_test(model, T = 300.0, erate = 1.0e-4, emax = 0.1,
    nsteps = 100):
  """
    Generate results for a tension test
  """
  res = drivers.uniaxial_test(model, erate, T = T, emax = emax, 
      nsteps = nsteps, full_results = True)
  
  data = np.concatenate((np.array(res['time'])[:,None], 
    np.array(res['temperature'])[:,None],
    np.array(res['strain']),
    np.array(res['stress'])), axis = 1)

  return data

def run_compare_test(model, target_times, target_temps, target_strains):
  """
    Run through a single test
  """
  driver = drivers.Driver_sd(model)
  for t,e,T in zip(target_times[1:], target_strains[1:], target_temps[1:]):
    driver.strain_step(e,t,T)

  data = np.concatenate((np.array(driver.t_int)[:,None], 
    np.array(driver.T_int)[:,None],
    np.array(driver.mechanical_strain_int),
    np.array(driver.stress_int)), axis = 1)

  return data

def rtt(tdir):
  """
    Load a model from tdir, load the time, temperature, and strain data
    from the csv file, run the test, and compare the stresses
  """
  model = parse.parse_xml(os.path.join(tdir, default_xml_file),
      default_model_name)
  input_data = np.loadtxt(os.path.join(tdir, default_result_name), delimiter = ',',
      skiprows = 1)

  output_data = run_compare_test(model, input_data[:,0], input_data[:,1], 
      input_data[:,2:8])

  return input_data[:,8:14], output_data[:,8:14]

def rtt_restart(tdir, xml_file = "xmodel.xml", model_name = "model"):
  """
    Load a model from tdir, save it to XML, reload it from XML,
    load the time, temperature, and strain data
    from the csv file, run the test, and compare the stresses
  """
  omodel = parse.parse_xml(os.path.join(tdir, default_xml_file),
      default_model_name)
  
  temp_dir = tempfile.mkdtemp()
  mfile = os.path.join(temp_dir,  xml_file)
  omodel.save(mfile, model_name)
  model = parse.parse_xml(mfile, model_name)

  input_data = np.loadtxt(os.path.join(tdir, default_result_name), delimiter = ',',
      skiprows = 1)

  output_data = run_compare_test(model, input_data[:,0], input_data[:,1], 
      input_data[:,2:8])

  return input_data[:,8:14], output_data[:,8:14]

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description = "Generate regression test data")
  parser.add_argument("directory", help = "Directory with the XML file.")
  parser.add_argument("--temperature", default = 300.0, type = float,
      help = "Temperature to run test at")
  parser.add_argument("--strain-rate", default = 1.0e-4, type = float,
      help = "Strain rate to run tension test")
  parser.add_argument("--maximum-strain", default = 0.1, type = float,
      help = "Maximum strain to target")
  parser.add_argument("--nsteps", default = 100, type = int,
      help = "Number of load steps to use")

  args = parser.parse_args()

  if not os.path.exists(args.directory):
    raise ValueError("Directory %s does not exist!" % args.directory)

  model_file = os.path.join(args.directory, default_xml_file)
  if not os.path.exists(model_file):
    raise ValueError("Model file %s does not exist!" % model_file)
  
  try:
    model = parse.parse_xml(model_file, default_model_name)
  except Exception:
    raise ValueError("Cannot load model name %s from file %s!" % (model_file, 
      default_model_name))

  result = generate_tensile_test(model, T = args.temperature, 
      erate = args.strain_rate, emax = args.maximum_strain,
      nsteps = args.nsteps)
  
  outfile = os.path.join(args.directory, default_result_name)
  np.savetxt(outfile, result, delimiter = ',', comments = '',
      header = "time,temperature,strain_11,strain_22,strain_33,strain_23,strain_13,strain_12,stress_11,stress_22,stress_33,stress_23,stress_13,stress_12")
