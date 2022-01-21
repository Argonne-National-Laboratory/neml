#include "pyhelp.h"

#include "math/rotations.h"

#include <iostream>
#include <string>
#include <stdexcept>

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml {

PYBIND11_MODULE(rotations, m) {
  py::module::import("neml.objects");

  m.doc() = "Orientation representations and related routines";

  py::class_<Quaternion, std::shared_ptr<Quaternion>>(m, "Quaternion")
      .def(py::init<const std::vector<double>>(), py::arg("vector"))
      .def("__repr__",
           [](Quaternion & me) -> std::string
           {
              std::ostringstream ss;
              ss << "Quaternion(array([" << me.quat()[0] << ", "
                << me.quat()[1] << ", " << me.quat()[2] << ", " <<
                me.quat()[3] << "]))";
              return ss.str();
           }, "python __repr__")

      .def("__str__",
           [](Quaternion & me) -> std::string
           {
              std::ostringstream ss;
              ss << "[" << me.quat()[0] << ", "
                << me.quat()[1] << ", " << me.quat()[2] << ", " <<
                me.quat()[3] << "]";
              return ss.str();
           }, "python __str__")
      .def_property_readonly("quat",
                             [](Quaternion & me) -> py::array_t<double>
                             {
                                auto vec = alloc_vec<double>(4);
                                std::copy(me.quat(), me.quat()+4, 
                                          arr2ptr<double>(vec));
                                return vec;
                             }, "Return the raw quaternion.")
      .def("norm", &Quaternion::norm, "2-norm of the components")
      .def("opposite", &Quaternion::opposite, "Opposite")
      .def("__neg__", &Quaternion::opposite, "Opposite with python sugar")
      .def("conj", &Quaternion::conj, "Conjugation")
      .def("flip", &Quaternion::flip, "Flip the scalar part")
      .def("inverse", &Quaternion::inverse, "Inverse")
      .def("exp", &Quaternion::exp, "Exponential map")
      .def("log", &Quaternion::log, "Logarithmic map")
      .def("pow", &Quaternion::pow, "Exponentiation")
      .def("__pow__", &Quaternion::pow, "Exponentiation with python sugar")
      .def("dot", &Quaternion::dot, "Cartesian dot product")
      .def_property_readonly("hash", &Quaternion::hash, "Hash for comparison")
      .def("to_product_matrix", 
           [](Quaternion & me) ->  py::array_t<double>
           {
            auto Mn = alloc_mat<double>(4,4);
            double * p = arr2ptr<double>(Mn);
            
            me.to_product_matrix(p);

            return Mn; 
            }, "Convert to a 4x4 matrix for making composition a dot product")

      .def(py::self *= py::self)
      .def(py::self *= double())
      .def(py::self * py::self)
      .def(double() * py::self)
      .def(py::self * double())
      
      .def(py::self /= py::self)
      .def(py::self /= double())
      .def(py::self / py::self)
      .def(py::self / double())
      ;

  py::class_<Orientation, Quaternion, std::shared_ptr<Orientation>>(m, "Orientation")
      .def(py::init<Quaternion&>(), "Copy from a quaternion object", py::arg("General quaternion"))
      .def(py::init(
           [](py::array_t<double> n, double a, std::string angles)
           {
            return new Orientation(Orientation::createAxisAngle(arr2ptr<double>(n), a,  angles));
           }), "Initialize from axis-angle representation",
           py::arg("axis"), py::arg("angle"), py::arg("angle_type") = "radians")
      .def(py::init(
           [](py::array_t<double> arr)
           {
            auto info = arr.request();
            if (info.ndim == 1) {
              if (info.shape[0] == 3) {
                return new Orientation(Orientation::createRodrigues(arr2ptr<double>(arr)));
              }
              else if (info.shape[0] == 4) {
                double * ptr = arr2ptr<double>(arr);
                return new Orientation(std::vector<double>(ptr,ptr+4));
              }
              else {
                throw std::runtime_error("Initializing a quaternion from a vector requires a length 3 or 4 array");
              }
            }
            else if (info.ndim == 2) {
              return new Orientation(Orientation::createMatrix(arr2ptr<double>(arr)));
            }
            else {
              throw std::runtime_error("Invalid buffer dimension, must be 1 or n");
            }
           }), "Initialize from a Rodrigues vector or a rotation matrix",
           py::arg("matrix"))
      .def(py::init(
           [](double a, double b, double c, std::string angles, std::string convention)
           {
            if (convention == "hopf") {
              return new Orientation(Orientation::createHopf(a, b, c, angles));
            }
            else if (convention == "hyperspherical") {
              return new Orientation(Orientation::createHyperspherical(a, b, c, angles));
            }
            else {
              return new Orientation(Orientation::createEulerAngles(a, b, c, angles, convention));
            }
           }), "Initialize from angles of various types",
           py::arg("a"), py::arg("b"), py::arg("c"), py::arg("angle_type") = "radians", py::arg("convention") = "kocks")
      .def(py::init(
              [](Vector & x, Vector & y)
              {
                return new Orientation(Orientation::createVectors(x,y));
              }), "Initialize from two orthogonal vectors.")
      .def("to_euler",
           [](Orientation & me, std::string angles, std::string convention) -> std::tuple<double, double, double>
           {
            double a, b, c;
            me.to_euler(a, b, c, angles, convention);
            return std::make_tuple(a,b,c);
           }, "Convert to Euler angles", py::arg("angle_type") = "radians", py::arg("convention") = "kocks")
      .def("to_axis_angle",
           [](Orientation & me, std::string angles) -> std::tuple<py::array_t<double>, double>
           {
            double n[3];
            double a;
            me.to_axis_angle(n, a, angles);
            auto nv = alloc_vec<double>(3);
            double * p = arr2ptr<double>(nv);
            std::copy(n, n+3, p);
            return std::make_tuple(nv, a);
           }, "Convert to axis-angle", py::arg("angle_type") = "radians")
      .def("to_matrix",
           [](Orientation & me) -> py::array_t<double>
           {
            double M[9];
            me.to_matrix(M);
            auto Mn = alloc_mat<double>(3,3);
            double * p = arr2ptr<double>(Mn);
            std::copy(M, M+9, p);
            return Mn;
           }, "Convert to a rotation matrix")
      .def("to_tensor", &Orientation::to_tensor)
      .def("to_rodrigues",
           [](Orientation & me) -> py::array_t<double>
           {
            double v[3];
            me.to_rodrigues(v);
            auto vn = alloc_vec<double>(3);
            double * p = arr2ptr<double>(vn);

            std::copy(v, v+3, p);

            return vn;
           }, "Convert to a Rodrigues vector")
      .def("to_hopf",
           [](Orientation & me, std::string angle) -> std::tuple<double, double, double>
           {
            double a, b, c;
            me.to_hopf(a, b, c, angle);
            return std::make_tuple(a,b,c);
           }, "Convert to Hopf coordinates", py::arg("angle_type") = "radians")
      .def("to_hyperspherical",
           [](Orientation & me, std::string angle) -> std::tuple<double, double, double>
           {
            double a1, a2, a3;
            me.to_hyperspherical(a1, a2, a3, angle);
            return std::make_tuple(a1, a2, a3);
           }, "Convert to hyperspherical coordinates", py::arg("angle_type") = "radians")

      .def("opposite", &Orientation::opposite, "Opposite")
      .def("__neg__", &Orientation::opposite, "Opposite with python sugar")
      .def("conj", &Orientation::conj, "Conjugation")
      .def("flip", &Orientation::flip, "Flip the scalar part")
      .def("inverse", &Orientation::inverse, "Inverse")
      .def("pow", &Orientation::pow, "Exponentiation")
      .def("__pow__", &Orientation::pow, "Exponentiation with python sugar")
      
      .def("apply", 
           [](Orientation & me, Vector & v) -> Vector
           {
            return me.apply(v);
           }, "Apply to a vector")
      .def("apply", 
           [](Orientation & me, RankTwo & v) -> RankTwo
           {
            return me.apply(v);
           }, "Apply to a RankTwo")
      .def("apply", 
           [](Orientation & me, Symmetric & v) -> Symmetric
           {
            return me.apply(v);
           }, "Apply to a Symmetric")
      .def("apply", 
           [](Orientation & me, Skew & v) -> Skew
           {
            return me.apply(v);
           }, "Apply to a Skew")
      .def("apply", 
           [](Orientation & me, RankFour & v) -> RankFour
           {
            return me.apply(v);
           }, "Apply to a RankFour")
      .def("apply", 
           [](Orientation & me, SymSymR4 & v) -> SymSymR4
           {
            return me.apply(v);
           }, "Apply to a SymSymR4")
      
      .def(py::self *= py::self)
      .def(py::self * py::self)
      
      .def(py::self /= py::self)
      .def(py::self / py::self)

      .def("distance", &Orientation::distance)
      ;

  py::class_<CrystalOrientation, Orientation, NEMLObject,
      std::shared_ptr<CrystalOrientation>>(m, "CrystalOrientation")
      .def(py::init(
           [](double a, double b, double c, std::string angles, std::string convention)
           {
            auto params = CrystalOrientation::parameters();
            params.assign_parameter("angles", std::vector<double>({a,b,c}));
            params.assign_parameter("angle_type", angles);
            params.assign_parameter("angle_convention", convention);

            return new CrystalOrientation(params);
           }), "Initialize from Euler angles",
           py::arg("a"), py::arg("b"), py::arg("c"), py::arg("angle_type") = "radians", py::arg("convention") = "kocks")
  ;
  
  m.def("random_orientations", &random_orientations);
  m.def("wexp", &wexp);
  m.def("wlog", &wlog);
  m.def("distance", &distance);
  m.def("rotate_to", &rotate_to);
  m.def("rotate_to_family", &rotate_to_family);
} // PYBIND11_MODULE(cpfmwk, m)

} // namespace cpfmwk
