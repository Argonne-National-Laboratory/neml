#!/usr/bin/env python

import numpy as np
import scipy.interpolate as inter
import scipy.optimize as opt
import os.path

import matplotlib.pyplot as plt

import sys
sys.path.append('../..')

from neml import drivers

path_to_xedl = '/home/messner/projects/xedl'
sys.path.append(path_to_xedl)

from query import merge_all
from xedl.xmlhelp import *

database = os.path.join(path_to_xedl,'database')

from generate import *

if __name__ == "__main__":
  # Setup the database
  data = merge_all(database)

  # Look at these temperatures
  tempsC = np.array([500, 550, 600])
  inc = 20.0

  for T in tempsC:
    Tk = T + 273.15
    # For each temperature find all the corresponding stress/strain
    # curves within inc degrees of the temperature
    experiments = data.xpath(
        'material[@name="gr91"]/experiment[@type="monotonic"][temperature<%f][temperature>%f][strain]'
        % (Tk + 20.0, Tk-20))
    # Summarize all the strain rates found in the experiments
    rates = sorted(list(set(float(e.find('rate').text) for e in experiments)))

    # Plot each experiment
    for experiment in experiments:
      strain = expand(experiment.find('strain').text)
      stress = expand(experiment.find('stress').text)
      plt.plot(strain, stress, 'k-')

    # Run and plot the Yaguchi model at each rate for this temperature
    yaguchi = make_yaguchi_model()
    for rate in rates:
      res = drivers.uniaxial_test(yaguchi, rate, T = Tk)
      plt.plot(res['strain'], res['stress'], 'r-')

    # Run and plot the Koo model at each rate for this temperature
    koo = make_koo_model(str(T))
    for rate in rates:
      res = drivers.uniaxial_test(koo, rate)
      plt.plot(res['strain'], res['stress'], 'b-')

    plt.show()

