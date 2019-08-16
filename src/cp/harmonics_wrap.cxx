#include "../pyhelp.h"

#include "harmonics.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml {

PYBIND11_MODULE(harmonics, m) {
  m.doc() = "Functions for computing spherical harmonics.";
  
  m.def("nharm_order", &nharm_order);
  m.def("nharm_total", &nharm_total);
  m.def("harmonic_SO3", &harmonic_SO3);
  m.def("P_SO3", &P_SO3);
  m.def("quadrature_S1", 
        [](size_t n) -> std::tuple<std::vector<double>, std::vector<double>>
        {
          std::vector<double> pts, wts;
          quadrature_S1(n, pts, wts);
          return std::make_tuple(pts, wts);
        }, "Quadrature over S1");
  m.def("quadrature_S2",
        [](size_t n) -> std::tuple<std::vector<double>, std::vector<double>,
        std::vector<double>>
        {
          std::vector<double> theta, phi, wts;
          quadrature_S2(n, theta, phi, wts);
          return std::make_tuple(theta, phi, wts);
        }, "Quadrature over S2");
  m.def("quadrature_SO3",
        [](size_t n) -> std::tuple<std::vector<Orientation>,
        std::vector<double>>
        {
          std::vector<Orientation> pts;
          std::vector<double> wts;
          quadrature_SO3(n, pts, wts);
          return std::make_tuple(pts, wts);
        }, "Quadrature over SO(3)");
  m.def("ds", &ds);
} // PYBIND11_MODULE

} // namespace neml
