#!/usr/bin/env python

from mpi4py import MPI

import numpy as np
import numpy.random as ra

from deap import base
from deap import creator
from deap import tools

import calibrate

class Member(list):
  def __init__(self, *args):
    list.__init__(self, *args)
    self.fitness = 0.0

def array2member(array):
  return [Member(list(a)) for a in array]

def member2array(pop):
  return np.array([np.array(a) for a in pop])

def assign_fitness(res, pop):
  for i, r in enumerate(res):
    pop[i].fitness = -r

  return pop

def log(f, ii, pop):
  residuals = [-m.fitness for m in pop]
  best = pop[np.argsort(residuals)[0]]

  f.write("Generation %i:\n" % (ii))
  f.write("Min:\t %e\n" % np.min(residuals))
  f.write("Max:\t %e\n" % np.max(residuals))
  f.write("Mean:\t %e\n" % np.mean(residuals))
  f.write("Best:\n")
  tw = ','.join(map(str, best))
  f.write(tw+'\n')
  f.write("\n")

def write_pop(fname, pop):
  res = [-m.fitness for m in pop]
  pop = np.array(pop)
  pop = pop[np.argsort(res)]
  np.save(fname, pop)

def fit_deap(model_maker, bounds, database, rweights, tweights, popsize = 50, 
    comm = MPI.COMM_WORLD, penalty = 10000.0, ngen = 500, cp = 0.5, 
    mp = 0.2, eta = 2.0, indpb = 0.05, ts = 3, elite = 5, 
    logfile = 'progress.log', popfile = 'population%i.npy',
    include_zero = False, population = None):
  """
    Use DEAP to fit a model to a database

    I use:
      two-point crossover with probability cp
      polynomial mutation with individual probability mp, attribute
        probability indpb, and parameter eta
      tournament selection with tournament size ts + elite elitism 

    Parameters:
      model_maker   function which takes a parameter list and returns a model
      bounds        bound on each parameter
      database      XEDL cases to fit to
      rweights      weights to apply to each residual
      tweights      takes a case and returns the approximate time

    Optional:
      popsize       population size
      comm          communicator to use
      penalty       how to penalize non-convergent results
      ngen          number of generations to run
      cp            crossover probability
      mp            mutation probability
      eta           mutation parameter
      indpb         independent probability of each gene being mutated
      ts            tournament size
      elite         elitism
      logfile       text file for loading information
      popfile       string formatter for generation files
      include_zero  include a special zero entry in the initial population
      population    specify the whole population (for restart)
  """
  rank = comm.Get_rank()
  
  # Some rearranging for DEAP
  ll = [l for (l,u) in bounds]
  ul = [u for (l,u) in bounds]

  # Root does most of this

  # Setup initial population
  if rank == 0:
    if population is not None:
      pop_array = population
    else:
      if include_zero:
        pop_array = np.array([[ra.uniform(p1,p2) for p1,p2 in bounds] for i in range(popsize-1)])
        pop_array = np.vstack((pop_array, np.zeros((len(bounds),))))
      else:
        pop_array = np.array([[ra.uniform(p1,p2) for p1,p2 in bounds] for i in range(popsize)])
    pop = array2member(pop_array)
    toolbox = base.Toolbox()
    toolbox.register("mate", tools.cxTwoPoint)
    toolbox.register("mutate", tools.mutPolynomialBounded, eta = eta, 
        low = ll, up = ul, indpb = indpb)
    toolbox.register("select", tools.selTournament, tournsize = ts)
    params = member2array(pop)
  else:
    params = None
 
  # evaluate initial residual
  res = evaluate_residual_mpi(model_maker, params, 
      database, rweights, tweights, comm = comm, penalty = penalty)
  #res = np.sum(np.abs(params), axis = 1)

  if rank == 0:
    # Assign
    pop = assign_fitness(res, pop)
    # Log
    log(open(logfile, 'w'), 0, pop)
    write_pop(popfile % 0, pop)


  # Loop over generations
  for ii in range(ngen):
    if rank == 0:
      # Start
      print("Generation %i" % (ii+1))

      # Select
      offspring = toolbox.select(pop, len(pop) - elite)
      offspring = list(map(toolbox.clone, offspring))

      # Crossover
      for child1, child2 in zip(offspring[::2], offspring[1::2]):
        if ra.random() < cp:
          toolbox.mate(child1, child2)

      # Mutate
      for mutant in offspring:
        if ra.random() < mp:
          toolbox.mutate(mutant)

      # Do our elitism
      iargs = np.argsort([m.fitness for m in pop])
      elites = [pop[i] for i in iargs[-elite:]]
      elites = list(map(toolbox.clone, elites))
      offspring += elites

      # Evaluate our population
      params = member2array(offspring)
      res = evaluate_residual_mpi(model_maker, params, database, rweights, tweights,
        comm = comm, penalty = penalty)
      #res = np.sum(np.abs(params), axis = 1)

      # Assign to the offspring
      offspring = assign_fitness(res, offspring)

      # Propagate the next generation
      pop[:] = offspring

      # Record some statistics
      residuals = [-o.fitness for o in pop]
      print("Best residual: %e" % np.min(residuals))
      print("")

      log(open(logfile, 'a'), ii+1, pop)
      write_pop(popfile % (ii+1), pop)


    else:
      params = None
      res = evaluate_residual_mpi(model_maker, params, database, rweights, tweights,
        comm = comm, penalty = penalty)

  return 

def evaluate_weights(database, tweights):
  """
    Evaluate the cost of each model run

    Parameters:
      database      database to fit to
      tweights      functions to determine cost
  """
  result = []
  for case in database:
    # Determine the type
    if case.attrib['type'] == 'monotonic' and case.find('control').text == 'strain':
      result.append(tweights['uniaxial'](case))
    elif case.attrib['type'] == 'relaxation' and case.find('control').text == 'strain':
      result.append(tweights['relax_strain'](case))
    elif case.attrib['type'] == 'relaxation' and case.find('control').text == 'stress':
      result.append(tweights['relax_stress'](case))
    elif case.attrib['type'] == 'cyclic' and case.find('control').text == 'strain':
      result.append(tweights['cyclic_strain'](case))
    elif case.attrib['type'] == 'cyclic' and case.find('control').text == 'stress':
      result.append(tweights['cyclic_stress'](case))
    else:
      raise ValueError("Unknown experiment with type %s." % case.attrib['type'])

  return result

def distribute_by_cost(costs, nmodels, size):
  """
    Do a round-robin distribution of the cases that need to be evaluated based
    on their cost

    Parameters:
      costs         cost of each evaluation
      nmodels       number of models we're going to run
      size          size of the communicator
  """
  distrib = -np.ones((nmodels,len(costs)), dtype = int)
  current = 0
  
  inds_cost = np.argsort(costs)

  for i in inds_cost:
    for j in range(nmodels):
      distrib[j,i] = current
      current += 1
      current %= size
  
  return distrib

def evaluate_residual_mpi(model_maker, params, database, rweights, tweights,
  comm = MPI.COMM_WORLD, penalty = 10000.0):
  """
    Actually evaluate the residual for a whole bunch of cases and models

    Parameters:
      model_maker   function which returns a model from parameters
      params        big list of parameters 
      database      list of XEDL databases
      rweights      weights to apply to each case
      tweights      approximate time for each evaluation

    Optional:
      comm          communicator to use
      penalty       penalty for "bad" results

    This function assumes that everyone has model_maker, database,
    and weights but only rank 0 has params

    Only root needs the final residual list
  """
  size = comm.Get_size()
  rank = comm.Get_rank()

  # Broadcast params as a numpy array
  if rank == 0:
    n = np.array(params.shape, dtype = int)
    comm.Bcast(n, root = 0)
    comm.Bcast(params, root = 0)
  else:
    n = np.zeros((2,), dtype = int)
    comm.Bcast(n, root = 0)
    params = np.zeros(n)
    comm.Bcast(params, root = 0)

  # Root runs weight calculations and does partitioning
  # Partitioning matrix is nmodels x ndatabase
  if rank == 0:
    costs = evaluate_weights(database, tweights)
    distrib = distribute_by_cost(costs, len(params), size)
  else:
    distrib = np.zeros((len(params), len(database)), dtype = int)
  
  # Broadcast the distribution matrix as a numpy array
  comm.Bcast(distrib, root = 0)

  # Now each proc can run it's cases.
  # Each proc must multiply the result by the appropriate weight before 
  # storing in the results array
  results = np.zeros((len(params), len(database)))

  # Not ideal, but hopefully this is relatively cheap
  models = [model_maker(p) for p in params]
  
  # Actually evaluate
  for i in range(distrib.shape[0]):
    for j in range(distrib.shape[1]):
      me = distrib[i,j]
      if me == rank:
        model = models[i]
        case = database[j]
        try:
          res = calibrate.evaluate_case(case, model, weights = rweights)
        except Exception as e:
          res = penalty
        if np.isnan(res) or np.isinf(res):
          res = penalty
        results[i,j] = res

  # Finally reduce to root and sum rows
  if rank == 0:
    rdata = np.zeros((len(params), len(database)))
    comm.Reduce(results, rdata, op = MPI.SUM, root = 0)
    return np.sum(rdata, axis = 1)
  else:
    comm.Reduce(results, None, op = MPI.SUM, root = 0)
    return None
