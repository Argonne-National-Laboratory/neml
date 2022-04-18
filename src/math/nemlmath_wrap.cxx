#include "pyhelp.h" // include first to avoid annoying redef warning

#include "math/nemlmath.h"

#include "nemlerror.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml {

PYBIND11_MODULE(nemlmath, m) {
  m.doc() = "Various mathematical helper functions.";
  
  m.def("transform_fourth",
        [](py::array_t<double> D, py::array_t<double> W) -> py::array_t<double>
        {
          if ((D.request().ndim != 2) || (D.request().shape[0] != 6) ||
              (D.request().shape[1] != 6)) {
            throw LinalgError("Input D is not a 6x6!");
          }
          if ((W.request().ndim != 2) || (W.request().shape[0] != 6) ||
              (W.request().shape[1] != 3)) {
            throw LinalgError("Input W is not a 6x3!");
          }

          auto A = alloc_mat<double>(9,9);

          transform_fourth(arr2ptr<double>(D), arr2ptr<double>(W), 
                           arr2ptr<double>(A));

          return A;

        }, "Transform the symmetric and skew parts to a full tensor");
  m.def("idsym",
        []()->py::array_t<double>
        {
          auto A = alloc_mat<double>(9,9);
          auto ptr = arr2ptr<double>(A);
          std::copy(idsym,idsym+81,ptr);
          return A;
        }, "The symmetric identity as a 9x9");
  m.def("idskew",
        []()->py::array_t<double>
        {
          auto A = alloc_mat<double>(9,9);
          auto ptr = arr2ptr<double>(A);
          std::copy(idskew,idskew+81,ptr);
          return A;
        }, "The skew symmetric identity as a 9x9");
  
  m.def("truesdell_tangent_outer",
        [](py::array_t<double, py::array::c_style> S) -> py::array_t<double>
        {
          if ((S.request().ndim != 1) || (S.request().shape[0] != 6)) {
            throw LinalgError("Input must be a Mandel vector");
          }

          auto A = alloc_mat<double>(9,9);
          truesdell_tangent_outer(arr2ptr<double>(S), arr2ptr<double>(A));

          return A;

        }, "Perform the outer product used in the Truesdell derivative and store as a 9x9");

  m.def("full2skew",
        [](py::array_t<double, py::array::c_style> A) -> py::array_t<double>
        {
          if (A.request().ndim != 2) {
            throw LinalgError("Input must be a 9x9 matrix!");
          }
          if ((A.request().shape[0] != 9) || (A.request().shape[1] != 9)) {
            throw LinalgError("Input must be a 9x9 matrix!");
          }

          auto M = alloc_mat<double>(6,3);
          full2skew(arr2ptr<double>(A), arr2ptr<double>(M));

          return M;
        }, "Convert a 9x9 to a skew-stored");
  m.def("skew2full",
        [](py::array_t<double, py::array::c_style> M) -> py::array_t<double>
        {
          if (M.request().ndim != 2) {
            throw LinalgError("Input must be a 6x3 matrix!");
          }
          if ((M.request().shape[0] != 6) || (M.request().shape[1] != 3)) {
            throw LinalgError("Input must be a 6x3 matrix!");
          }

          auto A = alloc_mat<double>(9,9);
          skew2full(arr2ptr<double>(M), arr2ptr<double>(A));

          return A;
        }, "Convert a skew-stored tensor to a 9x9");

  m.def("full2wws",
        [](py::array_t<double, py::array::c_style> A) -> py::array_t<double>
        {
          if (A.request().ndim != 2) {
            throw LinalgError("Input must be a 9x9 matrix!");
          }
          if ((A.request().shape[0] != 9) or (A.request().shape[1] != 9)) {
            throw LinalgError("Input must be a 9x9 matrix!");
          }

          auto M = alloc_mat<double>(3,6);
          full2wws(arr2ptr<double>(A), arr2ptr<double>(M));

          return M;
        }, "Convert a 9x9 to a ws matrix");
  
  m.def("wws2full",
        [](py::array_t<double, py::array::c_style> M) -> py::array_t<double>
        {
          if ((M.request().ndim != 2) || (M.request().shape[0] != 3) 
              || (M.request().shape[1] != 6)) {
            throw LinalgError("Input must be a 3x6 matrix!");
          }

          auto A = alloc_mat<double>(9,9);
          wws2full(arr2ptr<double>(M), arr2ptr<double>(A));

          return A;
        }, "Convert a ws to a full 9x9");

  m.def("full2mandel",
        [](py::array_t<double, py::array::c_style> A) -> py::array_t<double>
        {
          if (A.request().ndim != 2) {
            throw LinalgError("Input must be a 9x9 matrix!");
          }
          if ((A.request().shape[0] != 9) || (A.request().shape[1] != 9)) {
            throw LinalgError("Input must be a 9x9 matrix!");
          }

          auto M = alloc_mat<double>(6,6);
          full2mandel(arr2ptr<double>(A), arr2ptr<double>(M));

          return M;
        }, "Convert a 9x9 to a Mandel-stored");
  m.def("mandel2full",
        [](py::array_t<double, py::array::c_style> M) -> py::array_t<double>
        {
          if (M.request().ndim != 2) {
            throw LinalgError("Input must be a 6x6 matrix!");
          }
          if ((M.request().shape[0] != 6) || (M.request().shape[1] != 6)) {
            throw LinalgError("Input must be a 6x6 matrix!");
          }

          auto A = alloc_mat<double>(9,9);
          mandel2full(arr2ptr<double>(M), arr2ptr<double>(A));

          return A;
        }, "Convert a Mandel-stored tensor to a 9x9");
  m.def("truesdell_update_sym",
        [](py::array_t<double, py::array::c_style> D, py::array_t<double, py::array::c_style> W, 
           py::array_t<double, py::array::c_style> Sn, py::array_t<double, py::array::c_style> So) -> py::array_t<double>
        {
          if ((D.request().ndim != 1) || (W.request().ndim != 1) || (Sn.request().ndim != 1) || (So.request().ndim != 1)) {
            throw LinalgError("All inputs must be vectors!");
          }
          if ((D.request().shape[0] != 6) || (Sn.request().shape[0] != 6) || (So.request().shape[0]) != 6) {
            throw LinalgError("D, Sn, and So must be Mandel vectors!");
          }
          if (W.request().shape[0] != 3) {
            throw LinalgError("W must be a skew vector");
          }
          
          auto c = alloc_vec<double>(6);

          truesdell_update_sym(arr2ptr<double>(D), arr2ptr<double>(W), arr2ptr<double>(Sn), arr2ptr<double>(So),
                        arr2ptr<double>(c));

          return c;

        }, "Form the right hand side of the Truesdell convected update");
  m.def("truesdell_mat",
        [](py::array_t<double, py::array::c_style> D, py::array_t<double, py::array::c_style> W) -> py::array_t<double>
        {
          if ((D.request().ndim != 1) || (W.request().ndim != 1)) {
            throw LinalgError("All inputs must be vectors!");
          }
          if ((D.request().shape[0] != 6)) {
            throw LinalgError("D must be a Mandel vector!");
          }
          if (W.request().shape[0] != 3) {
            throw LinalgError("W must be a skew vector");
          }

          auto M = alloc_mat<double>(9,9);

          truesdell_mat(arr2ptr<double>(D), arr2ptr<double>(W), arr2ptr<double>(M));

          return M;

        }, "Form the matrix used in the convected update");
  m.def("truesdell_rhs",
        [](py::array_t<double, py::array::c_style> D, py::array_t<double, py::array::c_style> W, 
           py::array_t<double, py::array::c_style> Sn, py::array_t<double, py::array::c_style> So) -> py::array_t<double>
        {
          if ((D.request().ndim != 1) || (W.request().ndim != 1) || (Sn.request().ndim != 1) || (So.request().ndim != 1)) {
            throw LinalgError("All inputs must be vectors!");
          }
          if ((D.request().shape[0] != 6) || (Sn.request().shape[0] != 6) || (So.request().shape[0]) != 6) {
            throw LinalgError("D, Sn, and So must be Mandel vectors!");
          }
          if (W.request().shape[0] != 3) {
            throw LinalgError("W must be a skew vector");
          }
          
          auto c = alloc_vec<double>(6);

          truesdell_rhs(arr2ptr<double>(D), arr2ptr<double>(W), arr2ptr<double>(Sn), arr2ptr<double>(So),
                        arr2ptr<double>(c));

          return c;

        }, "Form the right hand side of the Truesdell convected update");
  
  m.def("sym",
        [](py::array_t<double, py::array::c_style> M) -> py::array_t<double>
        {
          if (M.request().ndim != 2) {
            throw LinalgError("Input must be a mtrix!");
          }
          if ((M.request().shape[0] != 3) || (M.request().shape[1] != 3)) {
            throw LinalgError("Input must be size 3x3!");
          }
          auto v = alloc_vec<double>(6);

          sym(arr2ptr<double>(M), arr2ptr<double>(v));

          return v;

        }, "Convert a full tensor to a Mandel vector.");

  m.def("usym",
        [](py::array_t<double, py::array::c_style> v) -> py::array_t<double>
        {
          if (v.request().ndim != 1) {
            throw LinalgError("Input must be a vector!");
          }
          if (v.request().shape[0] != 6) {
            throw LinalgError("Input must be a length 6 vector!");
          }
          auto M = alloc_mat<double>(3, 3);

          usym(arr2ptr<double>(v), arr2ptr<double>(M));
          
          return M;
        }, "Convert a Mandel vector to a full tensor.");

  m.def("skew",
        [](py::array_t<double, py::array::c_style> M) -> py::array_t<double>
        {
          if (M.request().ndim != 2) {
            throw LinalgError("Input must be a mtrix!");
          }
          if ((M.request().shape[0] != 3) || (M.request().shape[1] != 3)) {
            throw LinalgError("Input must be size 3x3!");
          }
          auto v = alloc_vec<double>(3);

          skew(arr2ptr<double>(M), arr2ptr<double>(v));

          return v;

        }, "Convert a full tensor to a skew vector.");

  m.def("uskew",
        [](py::array_t<double, py::array::c_style> v) -> py::array_t<double>
        {
          if (v.request().ndim != 1) {
            throw LinalgError("Input must be a vector!");
          }
          if (v.request().shape[0] != 3) {
            throw LinalgError("Input must be a length 3 vector!");
          }
          auto M = alloc_mat<double>(3, 3);

          uskew(arr2ptr<double>(v), arr2ptr<double>(M));
          
          return M;
        }, "Convert a skew vector to a full tensor.");

  m.def("minus_vec",
        [](py::array_t<double, py::array::c_style> a) -> py::array_t<double>
        {
          if (a.request().ndim != 1) {
            throw LinalgError("Array must be a vector!");
          }
          minus_vec(arr2ptr<double>(a), a.request().shape[0]);

          return a;
        }, "Negate a vector IN PLACE.");

  m.def("add_vec",
        [](py::array_t<double, py::array::c_style> a, py::array_t<double, py::array::c_style> b) -> py::array_t<double>
        {
          if ((a.request().ndim != 1) || (b.request().ndim != 1)) {
            throw LinalgError("Both arrays must be vectors!");
          }

          if (a.request().shape[0] != b.request().shape[0]) {
            throw LinalgError("Both vectors must have the same length!");
          }

          auto c = alloc_vec<double>(a.request().shape[0]);

          add_vec(arr2ptr<double>(a), arr2ptr<double>(b), a.request().shape[0], arr2ptr<double>(c));

          return c;

        }, "Add two vectors.");

  m.def("sub_vec",
        [](py::array_t<double, py::array::c_style> a, py::array_t<double, py::array::c_style> b) -> py::array_t<double>
        {
          if ((a.request().ndim != 1) || (b.request().ndim != 1)) {
            throw LinalgError("Both arrays must be vectors!");
          }

          if (a.request().shape[0] != b.request().shape[0]) {
            throw LinalgError("Both vectors must have the same length!");
          }

          auto c = alloc_vec<double>(a.request().shape[0]);

          sub_vec(arr2ptr<double>(a), arr2ptr<double>(b), a.request().shape[0], arr2ptr<double>(c));

          return c;

        }, "Subtract two vectors.");

  m.def("dot_vec",
        [](py::array_t<double, py::array::c_style> a, py::array_t<double, py::array::c_style> b) -> double
        {
          if ((a.request().ndim != 1) || (b.request().ndim != 1)) {
            throw LinalgError("Both arrays must be vectors!");
          }

          if (a.request().shape[0] != b.request().shape[0]) {
            throw LinalgError("Both vectors must have the same length!");
          }

          return dot_vec(arr2ptr<double>(a), arr2ptr<double>(b), a.request().shape[0]);

        }, "Compute a dot product between two vectors.");

  m.def("norm2_vec",
        [](py::array_t<double, py::array::c_style> a) -> double
        {
          if (a.request().ndim != 1) {
            throw LinalgError("The array must be a vector!");
          }

          return norm2_vec(arr2ptr<double>(a), a.request().shape[0]);
        }, "Compute the two norm of a vector.");

  m.def("normalize_vec",
        [](py::array_t<double, py::array::c_style> a) -> py::array_t<double>
        {
          if (a.request().ndim != 1) {
            throw LinalgError("The array must be a vector!");
          }
          normalize_vec(arr2ptr<double>(a), a.request().shape[0]);

          return a;
        }, "Normalize a vector IN PLACE.");

  m.def("dev_vec",
        [](py::array_t<double, py::array::c_style> a) -> py::array_t<double>
        {
          if (a.request().ndim != 1) {
            throw LinalgError("The array must be a vector!");
          }
          if (a.request().shape[0] != 6) {
            throw LinalgError("a must be a 6-vector!");
          }
          
          dev_vec(arr2ptr<double>(a));

          return a;

        }, "Make a vector deviatoric, IN PLACE.");

  m.def("outer_vec",
        [](py::array_t<double, py::array::c_style> a, py::array_t<double, py::array::c_style> b) -> py::array_t<double>
        {
          if ((a.request().ndim != 1)  || (b.request().ndim != 1)) {
            throw LinalgError("The arrays must be vectors!");
          }
          
          auto C = alloc_mat<double>(a.request().shape[0], b.request().shape[0]);

          outer_vec(arr2ptr<double>(a), a.request().shape[0], arr2ptr<double>(b), b.request().shape[0], arr2ptr<double>(C));

          return C;
        }, "Outer product of a x b.");

  m.def("outer_update",
        [](py::array_t<double, py::array::c_style> a, py::array_t<double, py::array::c_style> b, py::array_t<double, py::array::c_style> C) -> py::array_t<double>
        {
          if ((a.request().ndim != 1)  || (b.request().ndim != 1)) {
            throw LinalgError("The arrays must be vectors!");
          }

          if (C.request().ndim != 2) {
            throw LinalgError("C must be a matrix!");
          }

          if ((C.request().shape[0] != a.request().shape[0]) || (C.request().shape[1] != b.request().shape[0])) {
            throw LinalgError("C is not conformable with axb!");
          }
          
          outer_update(arr2ptr<double>(a), a.request().shape[0], arr2ptr<double>(b), b.request().shape[0], arr2ptr<double>(C));

          return C;
        }, "Rank 2 update C += a x b.");

  m.def("outer_update_minus",
        [](py::array_t<double, py::array::c_style> a, py::array_t<double, py::array::c_style> b, py::array_t<double, py::array::c_style> C) -> py::array_t<double>
        {
          if ((a.request().ndim != 1)  || (b.request().ndim != 1)) {
            throw LinalgError("The arrays must be vectors!");
          }

          if (C.request().ndim != 2) {
            throw LinalgError("C must be a matrix!");
          }

          if ((C.request().shape[0] != a.request().shape[0]) || (C.request().shape[1] != b.request().shape[0])) {
            throw LinalgError("C is not conformable with axb!");
          }
          
          outer_update_minus(arr2ptr<double>(a), a.request().shape[0], arr2ptr<double>(b), b.request().shape[0], arr2ptr<double>(C));

          return C;
        }, "Rank 2 update C -= a x b.");

  m.def("mat_vec",
        [](py::array_t<double, py::array::c_style> A, py::array_t<double, py::array::c_style> b) -> py::array_t<double>
        {
          if (A.request().ndim != 2) {
            throw LinalgError("A must be a matrix!");
          }
          if (b.request().ndim != 1) {
            throw LinalgError("b must be a vector!");
          }
          if (A.request().shape[1] != b.request().shape[0]) {
            throw LinalgError("A and b are not conformable!");
          }

          auto c = alloc_vec<double>(A.request().shape[0]);
          
          mat_vec(arr2ptr<double>(A), A.request().shape[0], arr2ptr<double>(b), b.request().shape[0], arr2ptr<double>(c));

          return c;

        }, "Matrix-vector product c = A.b.");

  m.def("mat_vec_trans",
        [](py::array_t<double, py::array::c_style> A, py::array_t<double, py::array::c_style> b) -> py::array_t<double>
        {
          if (A.request().ndim != 2) {
            throw LinalgError("A must be a matrix!");
          }
          if (b.request().ndim != 1) {
            throw LinalgError("b must be a vector!");
          }
          if (A.request().shape[0] != b.request().shape[0]) {
            throw LinalgError("A and b are not conformable!");
          }

          auto c = alloc_vec<double>(A.request().shape[1]);
          
          mat_vec_trans(arr2ptr<double>(A), A.request().shape[1], arr2ptr<double>(b), b.request().shape[0], arr2ptr<double>(c));

          return c;

        }, "Matrix-vector product c = A.T.b.");

  m.def("mat_mat",
        [](py::array_t<double, py::array::c_style> A, py::array_t<double, py::array::c_style> B) -> py::array_t<double>
        {
          if ((A.request().ndim != 2) || (B.request().ndim != 2)) {
            throw LinalgError("A and B must be matrices!");
          }
          if (A.request().shape[1] != B.request().shape[0]) {
            throw LinalgError("A and B must be conformal!");
          }

          auto C = alloc_mat<double>(A.request().shape[0], B.request().shape[1]);
          
          mat_mat(A.request().shape[0], B.request().shape[1], A.request().shape[1],
                  arr2ptr<double>(A), arr2ptr<double>(B), arr2ptr<double>(C));

          return C;
        }, "Matrix-matrix product C = A.B.");

  m.def("mat_mat_ABT",
        [](py::array_t<double, py::array::c_style> A, py::array_t<double, py::array::c_style> B) -> py::array_t<double>
        {
          if ((A.request().ndim != 2) || (B.request().ndim != 2)) {
            throw LinalgError("A and B must be matrices!");
          }
          if (A.request().shape[1] != B.request().shape[1]) {
            throw LinalgError("A and B must be conformal!");
          }

          auto C = alloc_mat<double>(A.request().shape[0], B.request().shape[0]);
          
          mat_mat_ABT(A.request().shape[0], B.request().shape[0], A.request().shape[1],
                  arr2ptr<double>(A), arr2ptr<double>(B), arr2ptr<double>(C));

          return C;
        }, "Matrix-matrix product C = A.B.T");

  m.def("invert_mat", 
        [](py::array_t<double, py::array::c_style> A) -> py::array_t<double>
        {
          if (A.request().ndim != 2) {
            throw LinalgError("Array is not a matrix!");
          }
          if (A.request().shape[0] != A.request().shape[1]) {
            throw LinalgError("Matrix is not square!");
          }
          invert_mat(arr2ptr<double>(A), A.request().shape[0]);
          // Should check non-singular
          return A;
        }, "Invert a matrix IN PLACE.");

   m.def("solve_mat",
        [](py::array_t<double, py::array::c_style> A, py::array_t<double, py::array::c_style> b) -> py::array_t<double>
        {
          if (A.request().ndim != 2) {
            throw LinalgError("A is not a matrix!");
          }
          if (A.request().shape[0] != A.request().shape[1]) {
            throw LinalgError("A is not square!");
          }
          if (b.request().ndim != 1) {
            throw LinalgError("b is not a vector!");
          }
          if (A.request().shape[0] != b.request().shape[0]) {
            throw LinalgError("A and b are not conformable!");
          }

          solve_mat(arr2ptr<double>(A), A.request().shape[0], arr2ptr<double>(b));

          return b;
        }, "Solve Ax=b.");

   m.def("condition",
        [](py::array_t<double, py::array::c_style> A) -> double
        {
          if (A.request().ndim != 2) {
            throw LinalgError("A is not a matrix!");
          }
          if (A.request().shape[0] != A.request().shape[1]) {
            throw LinalgError("A is not square!");
          }

          return condition(arr2ptr<double>(A), A.request().shape[1]);
        }, "Calculate the approximate condition number of A.");

   m.def("polyval", &polyval,
         "Evaluate a polynomial at x, highest order term first.");

   m.def("poly_from_roots", &poly_from_roots, 
         "Setup a polynomial from roots");

   m.def("differentiate_poly", &differentiate_poly, 
         "Differentiate a polynomial", py::arg("poly"), py::arg("n") = 1);

   m.def("eigenvalues_sym",
         [](py::array_t<double, py::array::c_style> s) -> std::tuple<double, double, double>
         {
           double vals[3];

           eigenvalues_sym(arr2ptr<double>(s), vals);

           return std::make_tuple(vals[0],vals[1],vals[2]);
         }, "Eigenvalues of a symmetric matrix.");

   m.def("eigenvectors_sym",
         [](py::array_t<double, py::array::c_style> s) -> py::array_t<double>
         {
           auto V = alloc_mat<double>(3,3);
           
           eigenvectors_sym(arr2ptr<double>(s), arr2ptr<double>(V));

           return V;
         }, "Eigenvectors of a symmetric matrix.");

   m.def("I1",
         [](py::array_t<double, py::array::c_style> s) -> double
         {
           return I1(arr2ptr<double>(s));
         }, "First principal invariant.");

   m.def("I2",
         [](py::array_t<double, py::array::c_style> s) -> double
         {
          return I2(arr2ptr<double>(s));
         }, "Second principal invariant.");

   m.def("dgttrf",
         [](py::array_t<double, py::array::c_style> DL, py::array_t<double, py::array::c_style> D, py::array_t<double, py::array::c_style> DU) ->
         std::tuple<py::array_t<double>, py::array_t<double>, py::array_t<double>, py::array_t<double>, py::array_t<int>>
         {
          auto DL_r = DL.request();
          auto D_r = D.request();
          auto DU_r = DU.request();

          if ((DL_r.ndim != 1) or 
              (D_r.ndim != 1) or
              (DU_r.ndim != 1)) {
            throw LinalgError("Diagonals should all be vectors!");
          }

          int n = D_r.shape[0];
          if ((DL_r.shape[0] != n-1) or (DU_r.shape[0] != n-1)) {
            throw LinalgError("Diagonals do not have compatible sizes!");
          }
          auto DU2 = alloc_vec<double>(n-2);
          auto IPIV = alloc_vec<int>(n);
          
          int ier;
          double * DLp = arr2ptr<double>(DL);
          double * Dp = arr2ptr<double>(D);
          double * DUp = arr2ptr<double>(DU);
          double * DU2p = arr2ptr<double>(DU2);
          int * IPIVp = arr2ptr<int>(IPIV);
          
          Py_BEGIN_ALLOW_THREADS
          dgttrf_(n, DLp, Dp, DUp, DU2p,IPIVp, ier);
          Py_END_ALLOW_THREADS
          
          if (ier != 0) {
            throw LinalgError("Diagonal factorization failed!");
          }

          return std::make_tuple(DL, D, DU, DU2, IPIV);

         }, "Thin wrapper for DGTTRF");

   m.def("dgttrs",
         [](py::array_t<double, py::array::c_style> DL, py::array_t<double, py::array::c_style> D, py::array_t<double, py::array::c_style> DU, py::array_t<double, py::array::c_style> DU2, py::array_t<int, py::array::c_style> IPIV, py::array_t<double, py::array::c_style> b) ->
         py::array_t<double>
         {
          auto DL_r = DL.request();
          auto D_r = D.request();
          auto DU_r = DU.request();
          auto DU2_r = DU2.request();
          auto IPIV_r = IPIV.request();
          auto b_r = b.request();

          if ((DL_r.ndim != 1) or 
              (D_r.ndim != 1) or
              (DU_r.ndim != 1) or
              (DU2_r.ndim != 1) or
              (IPIV_r.ndim != 1) or
              (b_r.ndim != 1)) {
            throw LinalgError("Diagonals should all be vectors!");
          }

          int n = D_r.shape[0];
          if ((DL_r.shape[0] != n-1) or (DU_r.shape[0] != n-1) 
              or (DU2_r.shape[0] != n-2) or (IPIV_r.shape[0] != n)
              or (b_r.shape[0] != n)) {
            throw LinalgError("Diagonals do not have compatible sizes!");
          }
          int ier;

          double * DLp = arr2ptr<double>(DL);
          double * Dp = arr2ptr<double>(D);
          double * DUp = arr2ptr<double>(DU);
          double * DU2p = arr2ptr<double>(DU2);
          int * IPIVp = arr2ptr<int>(IPIV);
          double * bp = arr2ptr<double>(b);
          
          Py_BEGIN_ALLOW_THREADS
          dgttrs_("N", n, 1, DLp, Dp, DUp, DU2p, IPIVp, bp, n, ier);
          Py_END_ALLOW_THREADS

          if (ier != 0) {
            throw LinalgError("Diagonal factorization failed!");
          }
          
          return b;

         }, "Thin wrapper for DGTTRS");

  m.def("gcd", &gcd);
  m.def("common_gcd", &common_gcd);
  m.def("reduce_gcd", &reduce_gcd);
  m.def("convert_angle", &convert_angle);
  m.def("cast_angle", &cast_angle);
  m.def("isclose", &isclose);
  m.def("rotate_matrix",
        [](py::array_t<double, py::array::c_style> A, py::array_t<double,
           py::array::c_style> B) -> py::array_t<double>
        {
          int m = A.request().shape[0];
          int n = A.request().shape[1];

          if ((A.request().ndim != 2) || (B.request().ndim != 2)) {
            throw LinalgError("A and B must be matrices!");
          }
          if ((n != B.request().shape[0]) || (n != B.request().shape[1])) {
            throw LinalgError("A and B must be conformal!");
          }

          auto C = alloc_mat<double>(m,m);
          
          rotate_matrix(m, n, arr2ptr<double>(A), arr2ptr<double>(B), arr2ptr<double>(C));

          return C; 
        }, "A * B * A.T");
  m.def("fact", &fact);
  m.def("factorial", &factorial);
}

} // namespace neml
