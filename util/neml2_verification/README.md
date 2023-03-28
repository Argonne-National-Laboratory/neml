This utility runs a model for comparison to NEML2.  A complete 
test case includes:

1. A NEML XML file defining the model.
2. A plaintext test file with the following format:

```
# These are comment lines
# Header: with the following space-separated entries
# 1. NEML XML file name
# 2. NEML model name
# 3. NEML2 serialization file name
# 4. NEML2 model name
# 5. Either "with_temperature" or "no_temperature"
# 6. rtol for comparison
# 7. atol for comparison
viscoplastic.xml viscoplastic_linear_isotropic xx xx no_temperature 1e-5 1e-8
# A line with a test description that could be read by a test harness
Compare Perzyna viscoplasticity with linear isotropic hardening
# Test input as a space-separated list of time strain stress and 
# (optionally) temperature values
# time strain_xx strain_yy strain_zz strain_yz strain_xz strain_xy stress_xx stress_yy stress_zz stress_yz stress_xz stress_xy (temperature)
0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
1.0 0.001 -0.0005 -0.0005 0.0 0.0 0.0 100 0 0 0 0 0 0
...
```

3. Eventually, a serialized NEML2 model file.

This utility takes the path to a NEML XML file, the name of the model, 
the strain rate, the maximum strain, the number of steps, and (optionally)
an isothermal temperature, runs the model, and writes the test file
in the same directory as the XML file.
