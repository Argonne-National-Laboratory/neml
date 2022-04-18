#include "pyhelp.h" // include first to avoid redef warning

#include "cp/batch.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml {

PYBIND11_MODULE(batch, m) {
  py::module::import("neml.cp.singlecrystal");
  py::module::import("neml.models");

  m.doc() = "Parallel batch evaluator for crystal models";

  m.def("evaluate_crystal_batch",
        [](SingleCrystalModel & model, 
           py::array_t<double, py::array::c_style> d_np1,
           py::array_t<double, py::array::c_style> d_n,
           py::array_t<double, py::array::c_style> w_np1,
           py::array_t<double, py::array::c_style> w_n,
           py::array_t<double, py::array::c_style> T_np1,
           py::array_t<double, py::array::c_style> T_n,
           double t_np1, double t_n,
           py::array_t<double, py::array::c_style> s_n,
           py::array_t<double, py::array::c_style> h_n,
           py::array_t<double, py::array::c_style> u_n,
           py::array_t<double, py::array::c_style> p_n,
           int nthreads ) ->
        std::tuple<
          py::array_t<double>, py::array_t<double>,
          py::array_t<double>, py::array_t<double>,
          py::array_t<double>, py::array_t<double>>
        {
          int n = d_np1.request().shape[0];
          if ((d_n.request().shape[0] != n) ||
              (w_np1.request().shape[0] != n) || 
              (w_n.request().shape[0] != n) ||
              (T_np1.request().shape[0] != n) ||
              (T_n.request().shape[0] != n) ||
              (s_n.request().shape[0] != n) ||
              (h_n.request().shape[0] != n) ||
              (u_n.request().shape[0] != n) ||
              (p_n.request().shape[0] != n)) {
            throw std::runtime_error("Inputs do not have the same first dimension!");
          }

          if ((d_np1.request().ndim != 2) || (d_np1.request().shape[1] != 6))
          {
            throw std::runtime_error("d_np1 does not have the right shape");
          }

          if ((d_n.request().ndim != 2) || (d_n.request().shape[1] != 6))
          {
            throw std::runtime_error("d_n does not have the right shape");
          }

          if ((w_np1.request().ndim != 2) || (w_np1.request().shape[1] != 3))
          {
            throw std::runtime_error("w_np1 does not have the right shape");
          }

          if ((w_n.request().ndim != 2) || (w_n.request().shape[1] != 3))
          {
            throw std::runtime_error("w_n does not have the right shape");
          }

          if (T_np1.request().ndim != 1)
          {
            throw std::runtime_error("T_np1 does not have the right shape");
          }

          if (T_n.request().ndim != 1)
          {
            throw std::runtime_error("T_n does not have the right shape");
          }

          if ((s_n.request().ndim != 2) || (s_n.request().shape[1] != 6))
          {
            throw std::runtime_error("s_n does not have the right shape");
          }

          int nh = model.nstore();

          if ((h_n.request().ndim != 2) || (h_n.request().shape[1] != nh))
          {
            throw std::runtime_error("h_n does not have the right shape");
          }

          if (u_n.request().ndim != 1)
          {
            throw std::runtime_error("u_n does not have the right shape");
          }

          if (p_n.request().ndim != 1)
          {
            throw std::runtime_error("p_n does not have the right shape");
          }

          auto s_np1 = alloc_mat<double>(n,6);
          auto h_np1 = alloc_mat<double>(n,nh);
          auto A_np1 = alloc_3d<double>(n,6,6); 
          auto B_np1 = alloc_3d<double>(n,6,3);
          auto u_np1 = alloc_vec<double>(n);
          auto p_np1 = alloc_vec<double>(n);

          // You need to grab the pointers before you release the GIL
          double * d_np1_ptr = arr2ptr<double>(d_np1);
          double * d_n_ptr = arr2ptr<double>(d_n);
          double * w_np1_ptr = arr2ptr<double>(w_np1);
          double * w_n_ptr = arr2ptr<double>(w_n);
          double * T_np1_ptr = arr2ptr<double>(T_np1);
          double * T_n_ptr = arr2ptr<double>(T_n);
          double * s_np1_ptr = arr2ptr<double>(s_np1);
          double * s_n_ptr = arr2ptr<double>(s_n);
          double * h_np1_ptr = arr2ptr<double>(h_np1);
          double * h_n_ptr = arr2ptr<double>(h_n);
          double * A_np1_ptr = arr2ptr<double>(A_np1);
          double * B_np1_ptr = arr2ptr<double>(B_np1);
          double * u_np1_ptr = arr2ptr<double>(u_np1);
          double * u_n_ptr = arr2ptr<double>(u_n);
          double * p_np1_ptr = arr2ptr<double>(p_np1);
          double * p_n_ptr = arr2ptr<double>(p_n);

          // bye bye GIL
          {
            py::gil_scoped_release release;
            evaluate_crystal_batch(model, n,
                                         d_np1_ptr, d_n_ptr, w_np1_ptr, w_n_ptr,
                                         T_np1_ptr, T_n_ptr, t_np1, t_n,
                                         s_np1_ptr, s_n_ptr,
                                         h_np1_ptr, h_n_ptr,
                                         A_np1_ptr, B_np1_ptr, 
                                         u_np1_ptr, u_n_ptr,
                                         p_np1_ptr, p_n_ptr, nthreads);
          }

          return std::make_tuple(s_np1, h_np1, A_np1, B_np1, u_np1, p_np1);

        }, "Batch update for a crystal model", py::arg("model"), 
      py::arg("d_np1"), py::arg("d_n"), py::arg("w_np1"), py::arg("w_n"),
      py::arg("T_np1"), py::arg("T_n"), py::arg("t_np1"), py::arg("t_n"),
      py::arg("s_n"), py::arg("h_n"), py::arg("u_n"), py::arg("p_n"), 
      py::arg("nthreads") = 1
      );

  m.def("init_history_batch",
        [](SingleCrystalModel & model, size_t n) -> py::array_t<double>
        {
          int nh = model.nstore();
          auto h = alloc_mat<double>(n,nh);
          init_history_batch(model, n, arr2ptr<double>(h));
          return h;
        }, "Batch initialize history");
  m.def("set_orientation_passive_batch",
        [](SingleCrystalModel & model, py::array_t<double,
           py::array::c_style> h, std::vector<Orientation> qs)
        {
          if (h.request().ndim != 2) {
            throw std::runtime_error("History should have dim = 2");
          }
          size_t n = h.request().shape[0];

          if ((size_t) h.request().shape[1] != model.nstore()) {
            throw std::runtime_error("History input not compatible with model");
          }

          if (qs.size() != n) {
            throw std::runtime_error("History and input orientations are not compatible!");
          }

          set_orientation_passive_batch(model, n, 
                                                  arr2ptr<double>(h),
                                                  qs);
        }, "Set passive orientations from a big vector");
  m.def("get_orientation_passive_batch",
        [](SingleCrystalModel & model, py::array_t<double, py::array::c_style>
           h) -> std::vector<Orientation>
        {
          std::vector<Orientation> res;

          if (h.request().ndim != 2) {
            throw std::runtime_error("History should have dim = 2");
          }
          size_t n = h.request().shape[0];
          
          if ((size_t) h.request().shape[1] != model.nstore()) {
            throw std::runtime_error("History input not compatible with model");
          }

          get_orientation_passive_batch(model, n, 
                                                  arr2ptr<double>(h),
                                                  res);
          return res;
        }, "Get passive orientations from a polycrystal.");

} // PYBIND11_MODULE(batch, m)

} // namespace neml
