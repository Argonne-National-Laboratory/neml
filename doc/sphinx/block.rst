Block evaluation functions
==========================

Overview
--------

This module provides helper routines for evaluating batches of NEML models
sequentially.  These routines can take advantage of:

1. Shared memory parallelism to evaluate the stress updates on many cores.
2. Nicely-aligned memory.
3. Explicit block matrix forms for converting to and from Mandel notation

Module description
------------------

.. doxygenfunction:: neml::block_evaluate

This function takes input in full tensor notation, internally converts
to Mandel notation using highly-efficient block matrix multiplications,
runs the stress update, using OpenMP if enabled, and reconverts to full
tensors using more block matrix multiplications.
