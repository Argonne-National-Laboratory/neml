#!/usr/bin/env python3

import argparse, time, pickle
import os.path
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from neml import models, parse, drivers

xml = os.path.join(os.path.dirname(__file__),"reference.xml")
models = os.path.join(os.path.dirname(__file__),"models.txt")

strain_rate = 1.0e-4

align = 50

def run_model(model_name, xml_file, repeat):
  model = parse.parse_xml(xml_file, model_name)
  
  fstr = "    Running %s..." % model_name
  print(fstr, end = '')

  start = time.time()
  for i in range(repeat):
    res = drivers.uniaxial_test(model, strain_rate)
  end = time.time()

  avg = (end-start) / repeat
  
  pad = " " * (align - len(fstr))
  print(pad+"%f s (avg)" % avg)

  return res

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="Run regression tests.")
  parser.add_argument('--regenerate', action = 'store_true', 
      help = "Regenerate the reference data")
  parser.add_argument('--xml-file',
      default = xml, help = "XML file with the models")
  parser.add_argument('--model-file',
      default = models, help = "Text file with the list to run")
  parser.add_argument('--single-model',
      help = "If set only run a single model given by the argument")
  parser.add_argument('--repeat', default = 1, type = int, 
      help = "Repeat each test several times for timing information")

  args = parser.parse_args()

  run = [line.strip() for line in open(args.model_file)]
  
  if args.single_model:
    if args.single_model in run:
      run = [args.single_model]
    else:
      raise ValueError("Requested single model %s is not in the list"
          "of valid models!" % args.single_model)
  
  print("Running models:")
  results = {model: run_model(model, args.xml_file, args.repeat)
      for model in run}

  if args.regenerate:
    for k, v in results.items():
      fname = k + ".pickle"
      pickle.dump(v, open(fname, 'wb'))


