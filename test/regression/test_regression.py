from neml import models, parse, drivers

import os.path

import unittest, pickle
import numpy as np

from common import *

xml = localize("regression/reference.xml")
models = localize("regression/models.txt")

strain_rate = 1.0e-4

def test_all_regression():
  names = [line.strip() for line in open(models)]
  for model in names:
    # Too expensive
    if model == "linearcp":
      continue
    yield check_regression, model

def check_regression(model):
  nmodel = parse.parse_xml(xml, model)
  res = drivers.uniaxial_test(nmodel, strain_rate)
  reference = pickle.load(open(os.path.join('test/regression',
    model+'.pickle'), 'rb'))
  assert(np.allclose(res['strain'], reference['strain']))
  assert(np.allclose(res['stress'], reference['stress']))
  assert(np.allclose(res['energy_density'], reference['energy_density']))
  assert(np.allclose(res['plastic_work'], reference['plastic_work']))

