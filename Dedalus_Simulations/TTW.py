# Load modules:

import numpy as np
import dedalus.public as d3
import logging
logger = logging.getLogger(__name__)
import numpy.random as random
from scipy.special import erf
import sys

# Command line input for parameters:

if len(sys.argv)>1:
  savename = str(sys.argv[1]).replace('.', '_')
else:
  savename = 'snapshots'

if len(sys.argv)>2:
  Ro = float(sys.argv[2])
else:
  Ro = 1.0

if len(sys.argv)>3:
  E = float(sys.argv[3])
else:
  E = 0.1

if len(sys.argv)>4:
  Pr = float(sys.argv[4])
else:
  Pr = 1.0

# Define parameters:

Nx, Nz = 1024, 128               # number of gridpoints
Lx, Lz = 6, 1                    # domain size
T = 20                           # runtime 
nu = 1e-4                        # horizontal mixing for stability, approx: (Lx/Nx)**2 

B = 2/Lx                         # background gradient in x
saves = 100                      # number of saves (excluding initial save)
Δt = min(T/saves/5, Lx/Nx)       # min of [enough timesteps to accurately save] and [dx/(2*U) from CFL]
dealias = 3/2
dtype = np.float64
timestepper = d3.RK222

# Build grids:

coords = d3.CartesianCoordinates('x', 'z')
dist = d3.Distributor(coords, dtype=dtype)
xbasis = d3.RealFourier(coords['x'], size=Nx, bounds=(-Lx/2, Lx/2), dealias=dealias)
zbasis = d3.ChebyshevT(coords['z'], size=Nz, bounds=(-Lz/2, Lz/2), dealias=dealias)

# Define fields:

u = dist.Field(name='u', bases=(xbasis,zbasis))
v = dist.Field(name='v', bases=(xbasis,zbasis))
w = dist.Field(name='w', bases=(xbasis,zbasis))
p = dist.Field(name='p', bases=(xbasis,zbasis))
b = dist.Field(name='b', bases=(xbasis,zbasis))

tau_p  = dist.Field(name='tau_p')
tau_u1 = dist.Field(name='tau_u1', bases=(xbasis))
tau_u2 = dist.Field(name='tau_u2', bases=(xbasis))
tau_v1 = dist.Field(name='tau_v1', bases=(xbasis))
tau_v2 = dist.Field(name='tau_v2', bases=(xbasis))
tau_w1 = dist.Field(name='tau_w1', bases=(xbasis))
tau_w2 = dist.Field(name='tau_w2', bases=(xbasis))
tau_b1 = dist.Field(name='tau_b1', bases=(xbasis))
tau_b2 = dist.Field(name='tau_b2', bases=(xbasis))

# Define equations:

x, z = dist.local_grids(xbasis, zbasis)
ex, ez = coords.unit_vector_fields(dist)

dx = lambda A: d3.Differentiate(A,coords[0])
dz = lambda A: d3.Differentiate(A,coords[1])

lift_basis = zbasis.derivative_basis(1)
lift = lambda A: d3.Lift(A, lift_basis, -1)

uz = dz(u) + lift(tau_u1)
vz = dz(v) + lift(tau_v1)
bz = dz(b) + lift(tau_b1)

zf = dist.Field(name='zf', bases=(zbasis))
zf['g'] = z

problem = d3.IVP([p, b, u, v, w, tau_p, tau_b1, tau_b2, tau_u1, tau_u2, tau_v1, tau_v2, tau_w1, tau_w2], namespace=locals())

problem.add_equation("dt(u) - E*dz(uz)    - nu*dx(dx(u)) + dx(p) - v + lift(tau_u2) = - Ro*(u*dx(u) + w*uz) - B*zf")
problem.add_equation("dt(v) - E*dz(vz)    - nu*dx(dx(v))         + u + lift(tau_v2) = - Ro*(u*dx(v) + w*vz)       ")
problem.add_equation("                                     dz(p) - b + lift(tau_w2) =   0                         ")
problem.add_equation("dt(b) - E/Pr*dz(bz) - nu*dx(dx(b))    + Ro*B*u + lift(tau_b2) = - Ro*(u*dx(b) + w*bz)       ")
problem.add_equation("                                 dx(u) + dz(w) + lift(tau_w1) =   0                         ")

problem.add_equation(" dz(u)(z=Lz/2) = 0")
problem.add_equation("dz(u)(z=-Lz/2) = 0")
problem.add_equation(" dz(v)(z=Lz/2) = 0")
problem.add_equation("dz(v)(z=-Lz/2) = 0")
problem.add_equation(" dz(b)(z=Lz/2) = 0")
problem.add_equation("dz(b)(z=-Lz/2) = 0")
problem.add_equation("     w(z=Lz/2) = 0")
problem.add_equation("     w(z=-Lz/2) + tau_p = 0")
problem.add_equation("      integ(p) = 0")

# Build solver and set stop time:

solver = problem.build_solver(timestepper)
solver.stop_sim_time = T*(1+1/(5*saves))

# Set IC:

def ss(z):
    return np.sin(z/np.sqrt(2))*np.sinh(z/np.sqrt(2))
def sc(z):
    return np.cos(z/np.sqrt(2))*np.sinh(z/np.sqrt(2))
def cs(z):
    return np.sin(z/np.sqrt(2))*np.cosh(z/np.sqrt(2))
def cc(z):
    return np.cos(z/np.sqrt(2))*np.cosh(z/np.sqrt(2))
def Det(z):
    return cc(z)**2+ss(z)**2

def K(z,z0):
    return -z+((cc(z0)+ss(z0))*cs(z)+(cc(z0)-ss(z0))*sc(z))/(np.sqrt(2)*Det(z0))
def K_p(z,z0):
    return -1+(cc(z0)*cc(z)+ss(z0)*ss(z))*(cc(z)-ss(z))/Det(z0)
def K_pp(z,z0):
    return ((cc(z0)+ss(z0))*sc(z)-(cc(z0)-ss(z0))*cs(z))/(np.sqrt(2)*Det(z0))

def b0(x):
    return erf(np.sqrt(np.pi)*x/2)
def db0dx(x):
    return np.exp(-x**2*np.pi/4)
def d2b0dx2(x):
    return -np.pi*x/2*np.exp(-x**2*np.pi/4)

z0 = 1/np.sqrt(4*E)

b['g'] = b0(x) - B*x
p['g'] = z * (b0(x) - B*x)
u['g'] = -np.sqrt(E) * K_pp(z/np.sqrt(E), z0) * db0dx(x)
v['g'] = -np.sqrt(E) * K(z/np.sqrt(E), z0) * db0dx(x)
w['g'] = E * K_p(z/np.sqrt(E), z0) * d2b0dx2(x)

# Define saves:

snapshots = solver.evaluator.add_file_handler(savename, sim_dt=T/saves)
snapshots.add_task(u, name='u', scales=1)
snapshots.add_task(v, name='v', scales=1)
snapshots.add_task(w, name='w', scales=1)
snapshots.add_task(b, name='b', scales=1)
snapshots.add_task(p, name='p', scales=1)
snapshots.add_task(dx(b) + B, name='M2', scales=1)

# set CFL condition:

CFL = d3.CFL(solver, initial_dt=Δt, cadence=1, safety=0.5, threshold=0.05, max_change=4, min_change=0.25, max_dt=T/saves/5, min_dt=1e-4)
CFL.add_velocity(Ro*(u*ex + w*ez))

# run simulation:

try:
  while solver.proceed:
    timestep = CFL.compute_timestep()
    solver.step(timestep)
    if (solver.iteration-1) % 100 == 0:
      logger.info('i=%i, t=%e, dt=%e' %(solver.iteration, solver.sim_time, timestep))
except:
  logger.info(':(')
finally:
  solver.log_stats()