#!/usr/bin/env python3 

from neml import models, parse, drivers

import numpy as np

import shutil
import argparse
import os.path

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            prog = "NEML2 verification test generator",
            description = "Generates verification test files to cross "
            "reference against the new NEML2 implementation")
    parser.add_argument("file_path")
    parser.add_argument("model_name")
    parser.add_argument("strain_rate", type = float)
    parser.add_argument("max_strain", type = float, nargs = '+')
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
    
    if len(args.max_strain) == 1:
        res = drivers.uniaxial_test(model, args.strain_rate, T = T,
                emax = args.max_strain[0], nsteps = args.nsteps, 
                full_results = True)
        time = np.array(res['time'])
        strain = np.array(res['strain'])
        stress = np.array(res['stress'])
        temperature = np.array(res['temperature'])
    elif len(args.max_strain) == 6:
        pass
    else:
        raise ValueError("max_strain should either be a single value, "
                "indicating a uniaxial strain-controlled test, " 
                "or a list of 6 numbers for a fully strain-controlled "
                "test")
    
    
    if use_T:
        data = np.hstack((time[:,None],strain,stress,temperature[:,None]))
    else:
        data = np.hstack((time[:,None],strain,stress))
    
    with open(os.path.join(basepath, args.test_file), 'w') as f:
        f.write("# NEML2/NEML verification test\n")
        f.write("%s %s x x %s %e %e\n" % (mfile, args.model_name, T_option, args.rtol, args.atol))
        f.write("# Comment line\n")
        f.write("\n")
        f.write("# The following uses the Mandel convention to store the strain and stress\n")
        data_header = "# time strain_xx strain_yy strain_zz strain_yz strain_xz strain_xy stress_xx stress_yy stress_zz stress_yz stress_xz stress_xy"
        if use_T:
            data_header += " temperature"
        f.write(data_header + "\n")
        np.savetxt(f, data)
