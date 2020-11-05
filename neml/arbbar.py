import numpy as np
import numpy.linalg as la
import networkx as nx

import warnings

import scipy.optimize as opt

from neml import uniaxial

class BarModel(nx.MultiGraph):
  """
    Model of an arbitrary bar network.

    Graph structure:
      I am an networkx graph

      Bars are edge data
      Nodes are rigid links

      DoFs are nodal displacements
      Equations are nodal equilibrium

      The specific input structure is a networkx graph.
      You may label the nodes however you want.
      Each edge must be a bar object, as defined below.

    Boundary conditions:
      Boundary conditions are at nodes
      They can be imposed forces or imposed displacements
      They are specified as functions of time

      Additionally, temperature histories can be assigned to 
      each bar using the bar object

    Loading:
      Calling solve(dt) advances the model from t_n to t_n + dt

    Solution history:
      Each bar maintains its own history internally, 
        as defined in the bar class
      Each node has properties:
          "displacements"
          "forces"

        that are maintained as the structure evolves

        The user does not need to provide these in the networkx model --
        they will be initialized in this object's __init__
      Additionally, the global properties:
          "time"
          "energy"
          "dissipation"
        are tracked 
  """
  def __init__(self, *args, **kwargs):
    super(nx.MultiGraph, self).__init__(*args, **kwargs)
    self.time = np.array([0.0])
    self.energy = np.array([0.0])
    self.dissipation = np.array([0.0])
    self.validated = False

    # Sign convention
    self.sign = [-1.0,1.0]

  def solve(self, dt, ndiv = 4, dfact = 2, verbose = False,
      extrapolate = False, rtol = 1.0e-6, atol = 1.0e-8):
    """
      Solve the next time increment

      Parameters:
        dt      time increment

      Optional:
        ndiv    max number of adaptive subdivisions
        dfact   what to divided by (SHOULD BE AN INTEGER)
        verbose dump a lot of solver info
    """
    if not self.validated:
      self.validate()
    
    free_nodes, fixed_nodes = self.free_fixed_nodes()

    dt_total = dt
    dt_attempt = dt
    step_total = dfact**ndiv
    step_attempt = dfact**ndiv
    step_curr = 0
    dtried = 0
    ef = 1.0

    while step_curr < step_total:
      solveme = lambda d: self.RJ(d, free_nodes, dt_attempt)
      d_guess = self.make_guess(free_nodes, extrapolate = extrapolate,
          value = ef)   
      try:
        d = newton(solveme, d_guess, verbose = verbose, fail_iter = True, 
            rtol = rtol, atol = atol)
        step_curr += step_attempt
      except Exception as e:
        dtried += 1
        if verbose:
          print("Substepping: %i" % dtried)
        if dtried == ndiv:
          raise e
        step_attempt /= dfact
        dt_attempt /= dfact
        ef /= dfact
        continue
      self.apply_update(d, free_nodes, dt_attempt)
      ef = 1.0

  def apply_update(self, d, free_nodes, dt):
    """
      Actually apply the results of an update

      Parameters:
        d           free displacements
        free_nodes  the list of free nodes
        dt          time increment
    """
    t_next = self.time[-1] + dt
    
    for n in self.nodes():
      # This gets summed below
      self.nodes[n]['forces'] = np.append(self.nodes[n]['forces'], 0.0)
      if n in free_nodes:
        self.nodes[n]['displacements'] = np.append(
            self.nodes[n]['displacements'], d[free_nodes.index(n)])
      else:
        self.nodes[n]['displacements'] = np.append(
            self.nodes[n]['displacements'], 
            self.nodes[n]['displacement bc'](t_next))
    
    e_sum = 0.0
    p_sum = 0.0
    for i,j,data in self.edges(data=True):
      me = data['object']
      for ni,n in enumerate((i,j)):
        self.nodes[n]['forces'][-1] += self.sign[ni] * me.stress_next * me.A
      e_sum += me.energy_next * me.A * me.l
      p_sum += me.dissipation_next * me.A * me.l
      me.advance_step()

    self.time = np.append(self.time, t_next)
    self.energy = np.append(self.energy, e_sum)
    self.dissipation = np.append(self.dissipation, e_sum)
  
  def make_guess(self, free_nodes, extrapolate = True, value = 1.0):
    """
      Make a displacement guess

      Parameters:
        free_nodes  the list of free nodes (this sets dof numbering)
    """
    if (not extrapolate) or (len(self.time) < 2):
      return np.array([self.nodes[n]['displacements'][-1] for n in free_nodes])
    else:
      xprev = np.array([self.nodes[n]['displacements'][-1] for n in free_nodes])
      xprevprev = np.array([self.nodes[n]['displacements'][-2] for n in free_nodes])
      return xprev + value * (xprev - xprevprev)
  
  def free_fixed_nodes(self):
    """
      Return lists of the free and fixed nodes
    """
    free = [n for n in self.nodes() if 'displacement bc' not in self.nodes[n]]
    fixed = [n for n in self.nodes() if 'displacement bc' in self.nodes[n]]

    return free, fixed

  def RJ(self, d, free_nodes, dt):
    """
      Compute the residual equilibrium equation and the Jacobian

      Parameters:
        d           current free displacements
        free_nodes  listing of free nodes (sets the dof order)
        dt          time increment
    """
    R = np.zeros((len(free_nodes),))
    J = np.zeros((len(free_nodes), len(free_nodes)))

    t_next = self.time[-1] + dt

    # External forces
    for i,n in enumerate(free_nodes):
      if 'force bc' in self.nodes[n]:
        R[i] -= self.nodes[n]['force bc'](t_next)
      
    for i,j,data in self.edges(data=True):
      # Delta d
      dd = 0.0
      for index, node in enumerate((i,j)):
        if node in free_nodes:
          dd += self.sign[index]*d[free_nodes.index(node)]
        else:
          dd += self.sign[index]*self.nodes[node]['displacement bc'](t_next)
      
      # Force update
      f, A = data['object'].update_force(dd, t_next, self.time[-1])

      # Assemble internal forces and derivative
      for index_i, node_i in enumerate((i,j)):
        if node_i in free_nodes:
          idx_i = free_nodes.index(node_i)
          R[idx_i] += self.sign[index_i] * f

          for index_j, node_j in enumerate((i,j)):
            if node_j in free_nodes:
              idx_j = free_nodes.index(node_j)
              J[idx_i, idx_j] += (self.sign[index_i]*self.sign[index_j] * A)
     
    return R, J

  def add_force_bc(self, node, bfn):
    """
      Add a force boundary condition

      Parameters:
        node        node to add to
        bfn         force as a function of time
    """
    if 'force bc' in self.nodes[node] or 'displacement bc' in self.nodes[node]:
      warnings.warn("Overriding previous boundary condition")
    if 'displacement bc' in self.nodes[node]:
      del self.nodes[node]['displacement bc']

    self.nodes[node]['force bc'] = bfn

  def add_displacement_bc(self, node, bfn):
    """
      Parameters:
        node        node to add to
        bfn         displacement as a function of time
    """
    if 'force bc' in self.nodes[node] or 'displacement bc' in self.nodes[node]:
      warnings.warn("Overriding previous boundary condition")
    if 'force bc' in self.nodes[node]:
      del self.nodes[node]['force bc']

    self.nodes[node]['displacement bc'] = bfn

  def validate(self):
    """
      Iterate through the graph and make sure that everything is okay to
      run.
    """
    for node in self.nodes():
      if 'displacements' not in self.nodes[node]:
        self.nodes[node]['displacements'] = np.array([0.0])
      if 'forces' not in self.nodes[node]:
        self.nodes[node]['forces'] = np.array([0.0])
    self.validated = True

  def gather_element(self, quantity):
    """
      Iterate through elements and gather a particular quantity

      Parameters:
        quantity        what to grab
    """
    return np.array([getattr(e['object'], quantity) for i,j,e in 
      self.edges(data = True)])

class Bar(object):
  """
    Object representing a single bar, including its history
  """
  def __init__(self, mat, A, l, T = lambda t: 0.0):
    """
      Parameters:
        mat     NEML small strain material model
        A       bar cross-sectional area
        l       bar length
      
      Optional:
        T       bar temperature as a function of time, defaults to 0 
    """
    self.mat = uniaxial.UniaxialModel(mat)
    self.A = A
    self.l = l
    self.T = T
    
    self.stress = np.zeros((1,))
    self.strain = np.zeros((1,))
    self.mstrain = np.zeros((1,))
    self.estrain = np.zeros((1,))
    self.tstrain = np.zeros((1,))
    self.energy = np.zeros((1,))
    self.dissipation = np.zeros((1,))
    self.temperature = np.array([self.T(0)])
    self.history = np.array([self.mat.init_store()])

    self.stress_next = 0.0
    self.strain_next = 0.0
    self.mstrain_next = 0.0
    self.estrain_next = 0.0
    self.tstrain_next = 0.0
    self.energy_next = 0.0
    self.dissipation_next = 0.0
    self.T_next = 0.0
    self.history_next = np.empty(self.history[0].shape)

  def update_force(self, d, t_np1, t_n):
    """
      Run a stress update, provide the force and the tangent

      Parameters:
        d       DELTA bar displacement
        t_np1   next time step
        t_n     last time step
    """
    self.T_next = self.T(t_np1)

    a_np1 = self.mat.alpha(self.T_next)
    a_n = self.mat.alpha(self.temperature[-1])
    dT = self.T_next - self.temperature[-1]

    self.tstrain_next = self.tstrain[-1] + dT * (a_np1 + a_n) / 2

    self.strain_next = d / self.l
    self.mstrain_next = self.strain_next - self.tstrain_next
    
    (self.stress_next, self.history_next, tangent, self.energy_next,
        self.dissipation_next) = self.mat.update(
            self.mstrain_next, self.mstrain[-1], 
            self.T_next, self.temperature[-1],
            t_np1, t_n,
            self.stress[-1], self.history[-1],
            self.energy[-1], self.dissipation[-1])

    self.estrain_next = self.mat.elastic_strains(self.stress_next, self.T_next,
        self.history_next)

    return self.stress_next * self.A, tangent * self.A / self.l

  def advance_step(self):
    """
      Actually append the speculative bar histories
    """
    self.stress = np.append(self.stress, self.stress_next)
    self.strain = np.append(self.strain, self.strain_next)
    self.mstrain = np.append(self.mstrain, self.mstrain_next)
    self.estrain = np.append(self.estrain, self.estrain_next)
    self.tstrain = np.append(self.tstrain, self.tstrain_next)
    self.energy = np.append(self.energy, self.energy_next)
    self.dissipation = np.append(self.dissipation, self.dissipation_next)
    self.temperature = np.append(self.temperature, self.T_next)
    self.history = np.vstack((self.history, self.history_next))
   
def newton(RJ, x0, verbose = False, rtol = 1.0e-6, atol = 1.0e-8, miter = 50,
    linesearch = 'none', bt_tau = 0.5, bt_c = 1.0e-4, fail_iter = True):
  """
    Manually-code newton-raphson so that I can output convergence info, if
    requested.

    Parameters:
      RJ        function return the residual + jacobian
      x0        initial guess

    Optional:
      verbose   verbose output
  """
  R, J = RJ(x0)
  nR = la.norm(R)
  nR0 = nR
  x = np.copy(x0)
  
  i = 0

  if verbose:
    print("Iter.\tnR\t\tnR/nR0\t\tcond\t\tlinesearch")
    print("%i\t%e\t%e\t" % (i, nR, nR / nR0))

  while (nR > rtol * nR0) and (nR > atol):
    a = la.solve(J, R)

    if linesearch == 'none':
      f = 1.0
    elif linesearch == 'backtracking':
      f = backtrack(RJ, R, J, x, -a, tau = bt_tau, c = bt_c, verbose = verbose)
    else:
      raise ValueError("Unknown linesearch type.")
    
    x -= (a * f)
    R, J = RJ(x)
    nR = la.norm(R)
    i += 1
    if verbose:
      print("%i\t%e\t%e\t%e\t%f" % (i, nR, nR / nR0,la.cond(J), f))
    if i > miter:
      if verbose:
        print("")
      if fail_iter:
        raise MaximumIterations()
      else:
        break
  
  if verbose:
    print("")

  return x

def backtrack(RJ, R0, J0, x, a, tau = 0.5, c = 1.0e-4, verbose = False):
  """
    Do a backtracking line search

    Parameters:
      RJ        residual/jacobian function
      R0        value of R, to save a feval
      J0        value of J, to save a feval
      x         point to start from
      a         direction

    Optional:
      tau       backtrack tau
      c         backtrack c
      verbose   verbose output
  """
  alpha = 1.0
  cv = la.norm(RJ(x + alpha * a)[0])
  while cv > la.norm(R0 + c * alpha * np.dot(J0, a)):
    alpha *= tau
    cv = la.norm(RJ(x + alpha * a)[0])

  return alpha
