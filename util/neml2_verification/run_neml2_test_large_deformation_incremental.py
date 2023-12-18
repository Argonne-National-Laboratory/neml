#!/usr/bin/env python3 

from neml import models, parse, drivers
from neml.math import nemlmath

import numpy as np
import numpy.linalg as la

import shutil
import argparse
import os.path

def drive_all_L(model, L, tmax, nsteps, T):
    """
        Strain controlled driver where the user gives a full 
        Mandel strain vector
    """
    d = nemlmath.sym(0.5*(L+L.T))
    w = nemlmath.skew(0.5*(L-L.T))
    times = np.linspace(0, tmax, nsteps+1)
    dt = tmax / nsteps

    h_n = model.init_store()
    
    stress = np.zeros((nsteps+1,6))
    d_n = np.zeros((6,))
    w_n = np.zeros((3,))
    temperature = np.zeros((nsteps+1,))
    temperature[0] = T

    u_n = 0.0
    p_n = 0.0
    t_n = 0.0


    for i in range(1,nsteps+1):
        d_np1 = d_n + d * dt
        w_np1 = w_n + w * dt
        t_np1 = t_n + dt
        temperature[i] = temperature[i-1]

        s_np1, h_np1, A_np1, B_np1, u_np1, p_np1 = model.update_ld_inc(d_np1, d_n, w_np1, w_n, temperature[i], temperature[i-1], t_np1, t_n, stress[i-1], h_n, u_n, p_n)

        d_n = np.copy(d_np1)
        w_n = np.copy(w_np1)
        stress[i] = s_np1
        h_n = np.copy(h_np1)
        t_n = t_np1
        u_n = u_np1
        p_n = p_np1

    deformation_rate = np.zeros((nsteps+1,6))
    deformation_rate[:] = d

    vorticity = np.zeros((nsteps+1,3))
    vorticity[:] = w

    return (times,
            deformation_rate, 
            vorticity,
            stress,
            temperature)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            prog = "NEML2 verification test generator",
            description = "Generates large deformation verification test files to cross "
            "reference against the new NEML2 implementation")
    parser.add_argument("file_path")
    parser.add_argument("model_name")
    parser.add_argument("spatial_velocity_gradient", type = float, nargs = 9)
    parser.add_argument("max_time", type = float)
    parser.add_argument("nsteps", type = int)
    parser.add_argument("test_file")
    parser.add_argument("--temperature", type = float,
            required = False)
    parser.add_argument("--rtol", type = float, 
            default = 1e-5)
    parser.add_argument("--atol", type = float,
            default = 1.0e-8)

    args = parser.parse_args()

    basepath = os.path.dirname(args.file_path)
    model = parse.parse_xml(args.file_path, args.model_name)
    mfile = os.path.basename(args.file_path)
    
    if args.temperature is not None:
        T = args.temperature
        use_T = True
        T_option = "with_temperature"
    else:
        T = 20.0
        use_T = False
        T_option = "no_temperature"
    
    # Here
    time, D, W, stress, temperature = drive_all_L(
            model, np.array(args.spatial_velocity_gradient).reshape(3,3), args.max_time, 
            args.nsteps, T)
    
    if use_T:
        data = np.hstack((time[:,None],D,W,stress,temperature[:,None]))
    else:
        data = np.hstack((time[:,None],D,W,stress))
    
    with open(os.path.join(basepath, args.test_file), 'w') as f:
        f.write("# NEML2/NEML verification test\n")
        f.write("%s %s x x %s %e %e\n" % (mfile, args.model_name, T_option, args.rtol, args.atol))
        f.write("# Comment line\n")
        f.write("\n")
        f.write("# The following uses the Mandel convention to store the strain and stress\n")
        data_header = "time deformation_rate_xx deformation_rate_yy deformation_rate_zz deformation_rate_yz deformation_rate_xz deformation_rate_xy vorticity_zy vorticity_xz vorticity_yx stress_xx stress_yy stress_zz stress_yz stress_xz stress_xy"
        if use_T:
            data_header += " temperature"
        f.write(data_header + "\n")
        np.savetxt(f, data)
