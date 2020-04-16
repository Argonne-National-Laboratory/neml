import numpy as np
import numpy.linalg as la
import scipy.linalg as sla
import scipy.linalg.lapack as lapack
from numpy.polynomial.legendre import leggauss as lgg
import scipy.optimize as opt

from neml import arbbar
from neml.math import nemlmath

from functools import partial

def generate_thickness_gradient(ri, ro, T1, T2, Tdot_hot, hold, 
    Tdot_cold = None, hold_together = 0.0, delay = 0.0):
  """
    Generate our usual through-wall temperature gradient

    Parameters:
      ri                inner radius
      ro                outer radius
      T1                initial temperature
      T2                hot temperature
      Tdot_hot          heating rate, also cooling rate if not specified
      hold              hold at gradient

    Optional:
      Tdot_cold         cooling rate
      hold_together     hold at no temperature gradient
  """
  if Tdot_cold is None:
    Tdot_cold = Tdot_hot

  dT = T2 - T1
  thot = np.abs(dT / Tdot_hot)
  tcold = np.abs(dT / Tdot_cold)
  period = thot + hold + tcold + hold_together

  def gradient(r, t):
    if t < delay:
      return T1
    t = (t-delay) % period

    if t < thot:
      Tright = T1 + t * Tdot_hot * np.sign(dT)
    elif t < thot + hold:
      Tright = T2
    elif t < thot + hold + tcold:
      Tright = T2 - (t - thot - hold) * Tdot_cold * np.sign(dT)
    else:
      Tright = T1

    xi = (r - ri) / (ro - ri)

    return (1-xi) * T1 + xi * Tright

  return gradient

def generate_standard_timesteps(T1, T2, Tdot_hot, hold, nheat, nhold,
    ncycles, 
    Tdot_cold = None, ncool = None, hold_together = 0.0,
    ntogether = None, delay = 0.0, ndelay = None):
  """
    Generate a standard time stepping scheme for the thickness gradients

    Parameters:
      T1                initial temperature
      T2                hot temperature
      Tdot_hot          heat rate, also cooling rate if not specified
      hold              hold at gradient
      nheat             number of steps during heating
      nhold             number of steps during hold
      ncycles           number of cycles

    Optional:
      Tdot_cold         cooling rate
      ncool             number of cooling steps
      hold_together     hold at zero gradient
      ntogether         number of iso hold steps
      delay             delay before cycling
      ndelay            number of delay steps
  """
  if Tdot_cold is None:
    Tdot_cold = Tdot_hot
  if ncool is None:
    ncool = nheat

  if (delay != 0.0) and (ndelay is None):
    raise ValueError("If the delay is not zero you need to provide ndelay!")

  if (hold_together != 0.0) and (ntogether is None):
    raise ValueError("If you provide hold_together you need to provide ntogether!")

  times = [0.0]
  if delay > 0.0:
    times.extend(list(times[-1] + np.linspace(0, delay, ndelay+1)[1:]))

  for i in range(ncycles):
    times.extend(list(times[-1] + np.linspace(0, np.abs(T2-T1) / Tdot_hot, 
      nheat+1)[1:]))
    if hold > 0.0:
      times.extend(list(times[-1] + np.linspace(0, hold, nhold+1)[1:]))
    times.extend(list(times[-1] + np.linspace(0, np.abs(T2-T1) / Tdot_cold,
      ncool+1)[1:]))
    if hold_together > 0.0:
      times.extend(list(times[-1] + np.linspace(0, hold_together, 
        ntogether+1)[1:]))

  return times

def generate_multimaterial_thickness_gradient(rs, ks, T1, T2, Tdot_hot,
    hold, Tdot_cold = None, hold_together = 0.0, delay = 0.0, flip = False):
  """
    Generate a thermal gradient alternating between hot and cold for the
    case of a multilayer material.

    This function assumes the "slab" limit, i.e. a thin walled section.

    Parameters:
      rs                radii
      ks                heat conduction coefficients
      T1                initial temperature
      T2                hot temperature
      Tdot_hot          heating rate, also cooling rate if not specified
      hold              hold at gradient

    Optional:
      Tdot_cold         cooling rate
      hold_together     hold at no temperature gradient
      flip              have gradient going the other direction
  """
  if Tdot_cold is None:
    Tdot_cold = Tdot_hot

  if len(rs) != (len(ks)+1):
    raise ValueError("Incompatible radii and thermal conductivity arrays")

  dT = T2 - T1
  thot = np.abs(dT / Tdot_hot)
  tcold = np.abs(dT / Tdot_cold)
  period = thot + hold + tcold + hold_together

  def gradient(r, t):
    if t < delay:
      return T1
    t = (t-delay) % period

    if t < thot:
      Tright = T1 + t * Tdot_hot * np.sign(dT)
    elif t < thot + hold:
      Tright = T2
    elif t < thot + hold + tcold:
      Tright = T2 - (t - thot - hold) * Tdot_cold * np.sign(dT)
    else:
      Tright = T1

    Rs = np.diff(rs) / np.array(ks)
    Rt = np.sum(Rs)
    
    if flip:
      Ti = list(T1 + (1-np.cumsum(np.insert(Rs,0,0.0)) / Rt) * (Tright - T1)) 
    else:
      Ti = list(T1 + np.cumsum(np.insert(Rs,0,0.0)) / Rt * (Tright - T1))

    return np.interp(r, rs, Ti)

  return gradient

def mesh(rs, ns, bias = False, n = 2.0):
  """
    Generate a 1D mesh

    Parameters
      rs        list of points which divide the line into segments
                each segment is guaranteed to be a node point
                
                x----------x-------x--------------x
                r0         r1      r2             r3

      ns        list giving the number of elements in each segment

    Optional:
      bias      bias the mesh towards the segment end points
      n         bias factor
  """
  if not bias:
    xpoints = [np.linspace(rs[0], rs[1], ns[0]+1)]
    for i in range(1, len(ns)):
      xpoints.append(np.linspace(rs[i], rs[i+1], ns[i]+1)[1:])
    xpoints = np.concatenate(tuple(xpoints))
  else:
    xpoints = [rs[0]]
    for i in range(len(rs)-1):
      if ns[i] % 2 == 0:
        even = True
        nuse = ns[i] // 2 + 1
      else:
        even = False
        nuse = (ns[i] + 1) // 2 + 1

      x1 = (np.linspace(0,1, nuse)[1:])**n
      x2 = 1-(np.linspace(0,1, nuse)[:-1])**n + 1.0
      if even:
        xss = np.concatenate((x1, np.flip(x2,0)))/2
      else:
        xss = np.concatenate((x1[:-1], np.flip(x2,0)))/2
      
      # Transform
      xss = xss * (rs[i+1] - rs[i]) + rs[i]
      xpoints.extend(list(xss))
    
    xpoints = np.array(xpoints)

  return xpoints

class VesselSectionProblem(object):
  """
    This serves as the superclass for both the models of a "section of a
    cylindrical vessel": the Bree problem and my steady-state axisymmetric 
    model.
  """
  def __init__(self, rs, mats, ns, T, p, rtol = 1.0e-6, atol = 1.0e-8,
      p_ext = lambda t: 0.0):
    """
      Parameters:
        rs      radii deliminating each region
        mats    material model for each region
        ns      number of elements per region
        T       temperature as a function of (r, t)
        p       pressure, as a function of (t)

      Optional:
        rtol    relative tolerance for N-R solve
        atol    absolute tolerance for N-R solve
        p_ext   external pressure, as a function of t
    """
    # Check
    if len(rs) - 1 != len(mats):
      raise ValueError("Length of radii must be one more then length of regions")
    if len(mats) != len(ns):
      raise ValueError("Inconsistent region definitions")

    self.rs = rs
    self.mats = mats
    self.ns = ns
    self.T = T
    self.p = p
    self.p_ext = p_ext

    self.rtol = rtol
    self.atol = atol

    self.ts = np.diff(rs)
    self.t = np.sum(self.ts)
    self.r = rs[-1]
    self.r_inner = rs[0]
    self.r_outer = rs[-1]

    self.times = [0.0]
    self.pressures = [self.p(0)]
    self.external_pressures = [self.p_ext(0.0)]

def rect(n):
  """
    Rectangle rule points and weights over the
    standard domain [-1,1]

    Parameters:
      n     number of points
  """
  npts = np.linspace(-1,1,n+1)
  pts = np.array([(a+b)/2 for a,b in zip(npts[0:],npts[1:])])
  weights = np.array([2.0/n]*len(pts))

  return pts, weights

class BreeProblem(VesselSectionProblem):
  """
    A Bree problem, referenced back to a vessel.

    If you'd rather have a "standard" Bree problem you can
    make the overall wall thickness = 1 spanning from
    r = [0,1] and directly apply the hoop stress as the 
    pressure.
  """
  def __init__(self, rs, mats, ns, T, p, rtol = 1.0e-6, atol = 1.0e-8,
      itype = "gauss", p_ext = lambda t: 0.0):
    """
      Parameters:
        rs      radii deliminating each region
        mats    material model for each region
        ns      number of elements per region
        T       temperature as a function of (r, t)
        p       pressure, as a function of (t)

      Optional:
        rtol    relative tolerance for N-R solve
        atol    absolute tolerance for N-R solve
        bias    bias the mesh towards the interfaces
        itype   "rectangle" or "gauss", defaults to "gauss"
        p_ext   external pressure as a function of t
    """
    super(BreeProblem, self).__init__(rs, mats, ns, T, p,
        rtol = rtol, atol = atol, p_ext = p_ext)
    
    self.regions = np.diff(self.rs)
    self.dlength = [r for r,n in zip(self.regions, self.ns) for i in range(n)]

    if itype == "gauss":
      ifn = lgg
    elif itype == "rectangle":
      ifn = rect
    else:
      raise ValueError("Unknown integration type %s!" % itype)

    self.ipoints = [(xi+1)*dl/2 + r for dl,r,n in zip(self.regions, self.rs, self.ns) for 
        xi in ifn(n)[0]]
    self.weights = [wi * dl/2 for dl,n in zip(self.regions, self.ns) for 
        wi in ifn(n)[1]]
    materials = [m for m,n in zip(self.mats, self.ns) for i in range(n)]

    self.npoints = len(self.ipoints)
    
    self.P = lambda t: (p(t) * self.r_inner - p_ext(t) * self.r_outer) / self.t

    # Setup the n bar model
    self.barmodel = arbbar.BarModel()
    self.barmodel.add_node(1)
    self.barmodel.add_node(2)

    for ipt, wt, mat in zip(self.ipoints, self.weights, materials):
      def tlocal(tt, x = ipt):
        return self.T(x, tt)
      self.barmodel.add_edge(1,2, object = arbbar.Bar(mat, wt, 1.0, 
        T = tlocal))

    self.barmodel.add_displacement_bc(1, lambda t: 0.0)
    self.barmodel.add_force_bc(2, lambda t: self.P(t) * np.sum(self.weights))

    # Setup fields
    self.axialstrain = [0.0]
    self.hoop = [self.P(0.0)]
    self.force = [0.0]
    self.temperatures = [np.array([self.barmodel[1][2][i]['object'].temperature[-1]
      for i in range(self.npoints)])]
    self.stresses = [np.array([self.barmodel[1][2][i]['object'].stress[-1]
      for i in range(self.npoints)])]
    self.tstrains = [np.array([self.barmodel[1][2][i]['object'].tstrain[-1]
      for i in range(self.npoints)])]
    self.mstrains = [np.array([self.barmodel[1][2][i]['object'].mstrain[-1]
      for i in range(self.npoints)])]
    self.estrains = [np.array([self.barmodel[1][2][i]['object'].estrain[-1]
      for i in range(self.npoints)])]
    self.histories = [np.array([self.barmodel[1][2][i]['object'].history[-1]
      for i in range(self.npoints)])]

    self.energy = [0.0]
    self.work = [0.0]

  def update_loading(self, new_p, new_T, new_p_ext = lambda t: 0.0):
    """
      Impose new loading functions
    """
    self.P = lambda t: (new_p(t) * self.r_inner - new_p_ext(t) * self.r_outer) / self.t
    self.barmodel.nodes[2]['force bc'] = lambda t: self.P(t) * np.sum(self.weights)

    self.T = new_T
    for i,ipt in enumerate(self.ipoints):
      def tlocal(tt, x = ipt):
        return self.T(x, tt) 
      self.barmodel[1][2][i]['object'].T = tlocal

  def step(self, t, rtol = 1.0e-6, atol = 1.0e-8, verbose = False,
      ndiv = 4, dfact = 2, extrapolate = False):
    """
      Advance to the next step

      Parameters:
        t           time requested

      Optional:
        rtol        solver relative tolerance
        atol        solver absolute tolerance
        verbose     print a lot of stuff
        ndiv        maximum number of adaptive subdivisions
        dfact       subdivision factor
        extrapolate try to extrapolate displacements
    """
    dt = t - self.times[-1]
    if dt < 0.0:
      raise ValueError("Requesting negative time step!")
    
    self.barmodel.solve(dt, ndiv = ndiv, dfact = dfact, 
        verbose = verbose, extrapolate = extrapolate)
    
    self.times.append(t)
    self.pressures.append(self.p(t))
    self.external_pressures.append(self.p_ext(t))
    self.hoop.append(self.P(t))
    self.force.append(self.barmodel.nodes[2]['forces'][-1])
    self.axialstrain.append(self.barmodel.nodes[2]['displacements'][-1])
    
    self.temperatures.append(np.array([self.barmodel[1][2][i]['object'].temperature[-1]
      for i in range(self.npoints)]))
    self.stresses.append(np.array([self.barmodel[1][2][i]['object'].stress[-1]
      for i in range(self.npoints)]))
    self.tstrains.append(np.array([self.barmodel[1][2][i]['object'].tstrain[-1]
      for i in range(self.npoints)]))
    self.mstrains.append(np.array([self.barmodel[1][2][i]['object'].mstrain[-1]
      for i in range(self.npoints)]))
    self.estrains.append(np.array([self.barmodel[1][2][i]['object'].estrain[-1]
      for i in range(self.npoints)]))
    self.histories.append(np.array([self.barmodel[1][2][i]['object'].history[-1]
      for i in range(self.npoints)]))

    self.energy.append(0.0)
    self.work.append(0.0)
    for i,wt in enumerate(self.weights):
      me = self.barmodel[1][2][i]['object']
      self.energy[-1] += wt*me.energy[-1]*me.A*me.l
      self.work[-1] += wt*me.dissipation[-1]*me.A*me.l
    
class AxisymmetricProblem(VesselSectionProblem):
  """
    Python driver for a axisymmetric reduced problem.

    My coordinate system uses:
    x = radius
    y = theta
    z = z

    when we call into the material model.

    The conditions here are "steady state axisymmetry."
    The model is 1D and represents the stress state away from structural 
    discontinuities.

    The default boundary conditions are axisymmetry + generalized plane strain.
    The generalized plane strain boundary condition prevents initially
    parallel cross sections from rotating, but allows some net axial strain.

    Using constrained = True gives traditional plane strain
  """
  def __init__(self, rs, mats, ns, T, p, rtol = 1.0e-6, atol = 1.0e-8,
      bias = False, factor = 2.0, ngpts = 1, constrained = False,
      p_ext = lambda t: 0.0):
    """
      Parameters:
        rs      radii deliminating each region
        mats    material model for each region
        ns      number of elements per region
        T       temperature as a function of (r, t)
        p       pressure, as a function of (t)

      Optional:
        rtol        relative tolerance for N-R solve
        atol        absolute tolerance for N-R solve
        bias        bias the mesh towards the interfaces
        factor      bias factor
        ngpts       number of gauss points per element
        constrained constrain the cylinder against thermal expansion
        p_ext       external pressure, as a function of time
        constrained fully constrain the cylinder against thermal expansion
    """
    super(AxisymmetricProblem, self).__init__(rs, mats, ns, T, p,
        rtol = rtol, atol = atol, p_ext = p_ext)

    self.constrained = constrained

    self.mesh = mesh(self.rs, self.ns, bias = bias, n = factor)
    self.nnodes = len(self.mesh)
    self.nelem = len(self.mesh) - 1
    self.ls = np.diff(self.mesh)

    self.materials = [m for m,n in zip(self.mats, self.ns) for i in range(n)]
    
    self.displacements = [[0.0 for i in range(self.nnodes)]]
    self.axialstrain = [0.0]
    
    self.ngpts = ngpts
    self.gpoints, self.gweights = lgg(ngpts)

    self.temperatures = [np.array([self.T(xi, 0) for xi in self.mesh])]

    self.histories = [[np.array([m.init_store() for xi in self.gpoints])
      for m in self.materials]]
    self.stresses = [np.array([np.zeros((self.ngpts,6)) for n in self.ns for i in range(n)])]
    self.strains = [np.array([np.zeros((self.ngpts,6)) for n in self.ns for i in range(n)])]
    self.tstrains = [np.array([np.zeros((self.ngpts,6)) for n in self.ns for i in range(n)])]
    self.mstrains = [np.array([np.zeros((self.ngpts,6)) for n in self.ns for i in range(n)])]
    self.estrains = [np.array([np.zeros((self.ngpts,6)) for n in self.ns for i in range(n)])]

    # Save a bit of time
    self.Nl = np.array([[(1-xi)/2, (1+xi)/2] for xi in self.gpoints])
    self.ri = np.array([np.dot(self.Nl, self.mesh[e:e+2]) for e in range(self.nelem)])
    self.Bl = np.array([[-1.0/2, 1.0/2] for xi in self.gpoints])

  def take_step(self, t, rtol = 1.0e-6, atol = 1.0e-8, ilimit = 10,
      verbose = False, predict = 1.0):
    """
      Actually take a time step, used to sort out adaptive integration
      
      Parameters:
        t       next time

      Optional:
        rtol    solver relative tolerance
        atol    solver absolute tolerance
        ilimit  solver iteration limit
        verbose flag to print a lot of convergence info
        predict factor for forward extrapolation for the initial guess:
                1.0 = full extrapolation
                0.0 = use previous step
    """
    T = np.array([self.T(xi, t) for xi in self.mesh])
    p_left = self.p(t)
    p_right = self.p_ext(t)

    # Get a guess
    if len(self.displacements) < 2:
      x0 = np.concatenate((self.displacements[-1] ,
        [self.axialstrain[-1]]))
    else:
      x0 = np.concatenate((self.displacements[-1] + (self.displacements[-1] - self.displacements[-2]) * predict,
        [self.axialstrain[-1] + (self.axialstrain[-1] - self.axialstrain[-2]) * predict]))

    x = np.copy(x0)
    if self.constrained:
      x[-1] = 0.0

    R, (J11_du,J11_d,J11_dl,J12,J21,J22), strains, tstrains, mstrains, estrains, stresses, histories = self.RJ(x, T, p_left,
        p_right, t)

    if self.constrained:
      nR0 = la.norm(R[:-1])
    else:
      nR0 = la.norm(R)
    nR = nR0

    if verbose:
      print("Iter.\tNorm res")
      print("%i\t%e" % (0, nR0))
    
    i = 0
    while (nR > atol) and (nR / nR0 > rtol) and (i < ilimit):
      if self.constrained:
        R1 = R[:-1]
        x1 = tri_solve(J11_du, J11_d, J11_dl, R1)
        x[:-1] -= x1

      else:
        R1 = R[:-1]
        R2 = R[-1]
        c = R1 - R2/J22 * J12
        u = -J12 / J22
        v = J21
        x1 = cute_solve(J11_du, J11_d, J11_dl, u, v, c)
        x2 = (R2 - np.dot(J21,x1)) / J22

        x[:-1] -= x1
        x[-1] -= x2

      R, (J11_du,J11_d,J11_dl,J12,J21,J22), strains, tstrains, mstrains, estrains, stresses, histories = self.RJ(x, T, p_left,
          p_right, t)

      if self.constrained:
        nR = la.norm(R[:-1])
      else:
        nR = la.norm(R)
      i += 1

      if verbose:
        print("%i\t%e" % (i, nR))

    if i == ilimit:
      raise RuntimeError("Exceeded maximum allowed iterations!")
    
    return p_left, p_right, T, strains, tstrains, mstrains, estrains, stresses, histories, x[:-1], x[-1]

  def step(self, t, rtol = 1.0e-6, atol = 1.0e-8, verbose = False,
      div = 2, max_div = 4, predict = 1.0, ilimit = 10):
    """
      Update to the next state

      Parameters:
        t       time requested

      Optional:
        rtol    solver relative tolerance
        atol    solver absolute tolerance
        verbose flag to print a lot of convergence info
        div     step division factor for adaptive stepping: 
                new_time_step = old_time_step / div
        max_div maximum number of adaptive subdivisions
        predict factor for forward extrapolation for the initial guess:
                1.0 = full extrapolation
                0.0 = use previous step
        ilimit  solver iteration limit
    """
    istep = div**max_div
    sstep = istep
    cstep = 0
    cdiv = 0
    dt = t - self.times[-1]

    ostep = len(self.times)

    while cstep != istep:
      cstep += sstep      
      try:
        p_left, p_right, T, strains, tstrains, mstrains, estrains, stresses, histories, disp, axial = self.take_step(
            self.times[-1] + dt * float(cstep)/istep, rtol = rtol, atol = atol,
            verbose = verbose, predict = predict, ilimit = ilimit)
        self.times.append(t)
        self.pressures.append(p_left)
        self.external_pressures.append(p_right)
        self.temperatures.append(T)
        self.strains.append(strains)
        self.tstrains.append(tstrains)
        self.mstrains.append(mstrains)
        self.estrains.append(estrains)
        self.stresses.append(stresses)
        self.histories.append(histories)
        self.displacements.append(disp)
        self.axialstrain.append(axial)
      except RuntimeError as e:
        cstep -= sstep
        cdiv += 1
        if cdiv > max_div:
          self.times = self.times[:ostep]
          self.pressures = self.pressures[:ostep]
          self.external_pressures = self.external_pressures[:ostep]
          self.temperatures = self.temperatures[:ostep]
          self.strains = self.strains[:ostep]
          self.tstrains = self.tstrains[:ostep]
          self.mstrains = self.mstrains[:ostep]
          self.estrains = self.estrains[:ostep]
          self.stresses = self.stresses[:ostep]
          self.histories = self.histories[:ostep]
          self.displacements = self.displacements[:ostep]
          self.axialstrain = self.axialstrain[:ostep]
          raise e
        sstep /= div
        if verbose:
          print("Substepping: new step fraction %f" % (float(sstep) / istep))

    self.times = self.times[:ostep] + [self.times[-1]]
    self.pressures = self.pressures[:ostep] + [self.pressures[-1]]
    self.external_pressures = self.external_pressures[:ostep] + [self.external_pressures[-1]]
    self.temperatures = self.temperatures[:ostep] + [self.temperatures[-1]]
    self.strains = self.strains[:ostep] + [self.strains[-1]]
    self.tstrains = self.tstrains[:ostep] + [self.tstrains[-1]]
    self.mstrains = self.mstrains[:ostep] + [self.mstrains[-1]]
    self.estrains = self.estrains[:ostep] + [self.estrains[-1]]
    self.stresses = self.stresses[:ostep] + [self.stresses[-1]]
    self.histories = self.histories[:ostep] + [self.histories[-1]]
    self.displacements = self.displacements[:ostep] + [self.displacements[-1]]
    self.axialstrain = self.axialstrain[:ostep] + [self.axialstrain[-1]]


  def strain(self, d, l, r, ez):
    """
      Compute the full strains

      Parameters:
        d   displacements
        l   length
        r   radius
        ez  axial strain
    """
    strain = np.zeros((len(r),6))
    strain[:,0] = np.dot(self.Bl, d) * 2.0 / l
    strain[:,1] = np.dot(self.Nl, d) / r
    strain[:,2] = ez 
    return strain

  def RJ(self, x, T, p_left, p_right, t):
    """
      Compute the residual, jacobian, and the rest of the updated quantities

      Parameters:
        x           iterate
        T           nodal temperatures
        Tn          previous temperatures
        p_left      left pressure
        t           time
        p_right     right pressure
        t           time 
    """
    d = x[:-1]
    ez = x[-1]

    Fext = np.zeros((self.nnodes+1,))
    Fext[0] = p_left
    Fext[-2] = -p_right
    #Fext[-1] = (p_left *self.r_inner - p_right * self.r_outer) / (2*self.t)
    # We should use the proper thick-walled vessel formula
    Fext[-1] = (p_left * self.r_inner**2.0 / (self.r_outer**2.0 - self.r_inner**2.0) 
        -p_right * self.r_outer**2.0 / (self.r_outer**2.0 - self.r_inner**2.0))
    
    strains = np.zeros((self.nelem, self.ngpts, 6))
    tstrains = np.zeros((self.nelem, self.ngpts, 6))
    mstrains = np.zeros((self.nelem, self.ngpts, 6))
    estrains = np.zeros((self.nelem, self.ngpts, 6))
    stresses = np.zeros((self.nelem, self.ngpts, 6))
    histories = []

    Fint = np.zeros((self.nnodes+1,))
    A11_dl = np.zeros((self.nnodes-1,))
    A11_d = np.zeros((self.nnodes,))
    A11_du = np.zeros((self.nnodes-1,))
    A12 = np.zeros((self.nnodes,))
    A21 = np.zeros((self.nnodes,))
    A22 = 0.0

    for e, (mat, l, rs) in enumerate(zip(self.materials, self.ls, self.ri)):
      de = d[e:e+2]
      T_e = np.dot(self.Nl, T[e:e+2])
      T_n = np.dot(self.Nl, self.temperatures[-1][e:e+2])
      
      therm_inc = np.array([mat.alpha(Ti)*(Ti-Tii)*np.array([1,1,1,0,0,0]) 
        for Ti, Tii in zip(T_e, T_n)])

      strain = self.strain(de, l, rs, ez)
      tstrain = self.tstrains[-1][e] + therm_inc
      mstrain = strain - tstrain
      mstrain_n = self.mstrains[-1][e]
      stress_n = self.stresses[-1][e]
      hist_n = self.histories[-1][e]

      stress = np.zeros((self.ngpts,6))
      history = np.zeros((self.ngpts,len(hist_n[0])))
      estrain = np.zeros((self.ngpts,6))

      for i in range(self.ngpts):
        si, hi, ti, ui, pi = mat.update_sd(mstrain[i], mstrain_n[i],
            T_e[i], T_n[i], t, self.times[-1], stress_n[i], hist_n[i], 
            0.0, 0.0)
        stress[i] = si
        history[i] = hi

        estrain[i] = mat.elastic_strains(si, T_e[i], hi)
        
        wi = self.gweights[i]
        #Actually add in our terms
        Fint[e:e+2] += wi * (si[0] * self.Bl[i] * 2.0 / l + (si[1] - si[0]) * self.Nl[i] / rs[i]) * l / 2.0
        Fint[-1] += wi * si[2] * l / 2.0 / self.t
        
        DE = np.array([self.Bl[i] * 2.0/l, self.Nl[i] / rs[i]])
        
        # Remember wi and the jacobian
        J11 = np.dot(np.outer(2.0/l * self.Bl[i], ti[0,:2]) + np.outer(self.Nl[i] / rs[i], 
            ti[1,:2] - ti[0,:2]), DE) * wi * l / 2.0 
        J12 = (ti[0,2] * self.Bl[i] * 2.0 / l + (ti[1,2] - ti[0,2]) * self.Nl[i] / rs[i]) * wi * l / 2.0
        J21 = (np.dot(ti[2,:2], DE) / self.t) * wi * l / 2.0 
        J22 = (ti[2,2] / self.t) * wi * l / 2.0

        A11_du[e] += J11[0,1]
        A11_d[e] += J11[0,0]
        A11_d[e+1] += J11[1,1]
        A11_dl[e] += J11[1,0]
        A12[e:e+2] += J12
        A21[e:e+2] += J21
        A22 += J22

      strains[e] = strain
      tstrains[e] = tstrain
      mstrains[e] = mstrain
      estrains[e] = estrain
      stresses[e] = stress
      histories.append(history)

    return Fint - Fext, (A11_du,A11_d,A11_dl,A12,A21,A22), strains, tstrains, mstrains, estrains, stresses, histories

def tri_solve(DU, D, DL, c):
  """
    Solve B x = c
  """
  DL, D, DU, DU2, IPIV = nemlmath.dgttrf(DL, D, DU)
  x = nemlmath.dgttrs(DL,D,DU,DU2,IPIV,c)

  return x

def cute_solve(DU, D, DL, u, v, c):
  """
    Do a cute solve of:

    (B + u v.T) x = c

    with B tri-diagonal
  """
  DL, D, DU, DU2, IPIV = nemlmath.dgttrf(DL, D, DU)
  
  x = nemlmath.dgttrs(DL,D,DU,DU2,IPIV, c)
  y = nemlmath.dgttrs(DL,D,DU,DU2,IPIV, u)
  
  return x - np.dot(v,x) / (1+np.dot(v,y)) * y
