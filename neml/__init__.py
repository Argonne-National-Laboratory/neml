import os

ndir = os.path.dirname(__file__)
os.environ['PATH'] += os.pathsep + ndir

if os.name == 'nt':
    os.add_dll_directory(ndir)

import neml.history
