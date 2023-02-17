import sys
sys.path.append('..')

from common import *

from neml import nemlerror
import unittest

class TestHandleException(unittest.TestCase):
    def test_NEMLError(self):
        self.assertRaises(nemlerror.NEMLError, 
                lambda: nemlerror.throw_exception(nemlerror.ExceptionType.NEMLError))

    def test_LinalgError(self):
        self.assertRaises(nemlerror.LinalgError, 
                lambda: nemlerror.throw_exception(nemlerror.ExceptionType.LinalgError))

    def test_NonlinearSolverError(self):
        self.assertRaises(nemlerror.NonlinearSolverError, 
                lambda: nemlerror.throw_exception(nemlerror.ExceptionType.NonlinearSolverError))
