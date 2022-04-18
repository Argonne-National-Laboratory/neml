#include "pyhelp.h" // avoid redef warning

#include "block.h"

#include "nemlerror.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml {

PYBIND11_MODULE(block, m) {
  py::module::import("neml.objects");

  m.doc() = "Functions for evaluating whole blocks of models at once";

  m.def("t2m_array",
        []() -> py::array_t<double>
        {
          auto array = alloc_tensor<double>({3,3,6});
          auto ptr = arr2ptr<double>(array);

          std::copy(t2m_array, t2m_array+54, ptr);

          return array;
        }, "Static tensor to Mandel array");
  m.def("m2t_array",
        []() -> py::array_t<double>
        {
          auto array = alloc_tensor<double>({6,3,3,});
          auto ptr = arr2ptr<double>(array);

          std::copy(m2t_array, m2t_array+54, ptr);

          return array;
        }, "Static Mandel to tensor array");
  m.def("m42t4_array",
        []() -> py::array_t<double>
        {
          auto array = alloc_tensor<double>({6,6,3,3,3,3});
          auto ptr = arr2ptr<double>(array);

          std::copy(m42t4_array, m42t4_array+2916, ptr);

          return array;
        }, "Static Mandel rank 4 to full rank 4 array");

  m.def("block_evaluate",
        [](std::shared_ptr<NEMLModel> model, 
           py::array_t<double, py::array::c_style> e_np1,
           py::array_t<double, py::array::c_style> e_n,
           py::array_t<double, py::array::c_style> T_np1,
           py::array_t<double, py::array::c_style> T_n,
           double t_np1, double t_n,
           py::array_t<double, py::array::c_style> s_np1,
           py::array_t<double, py::array::c_style> s_n,
           py::array_t<double, py::array::c_style> h_np1,
           py::array_t<double, py::array::c_style> h_n,
           py::array_t<double, py::array::c_style> A_np1,
           py::array_t<double, py::array::c_style> u_np1,
           py::array_t<double, py::array::c_style> u_n,
           py::array_t<double, py::array::c_style> p_np1,
           py::array_t<double, py::array::c_style> p_n
           )
        {
          auto gets = [](py::array_t<double> a, size_t d) -> size_t {return a.request().shape[d];};
          auto getn = [](py::array_t<double> a) -> size_t {return a.request().shape[0];};
          auto getd = [](py::array_t<double> a) -> size_t {return
            a.request().ndim;};
          size_t nblock = getn(e_np1);
          size_t nhist = model->nstore();

          // Check that all the arrays have the right leading dimension
          if ((getn(e_np1) != nblock) ||
              (getn(e_n) != nblock) ||
              (getn(T_np1) != nblock) ||
              (getn(T_n) != nblock) ||
              (getn(s_np1) != nblock) ||
              (getn(s_n) != nblock) ||
              (getn(h_np1) != nblock) ||
              (getn(h_n) != nblock) ||
              (getn(A_np1) != nblock) ||
              (getn(u_np1) != nblock) ||
              (getn(u_n) != nblock) ||
              (getn(p_np1) != nblock) ||
              (getn(p_n) != nblock)) {
            throw std::runtime_error("Inputs do not all have the same leading"
                                     "dimension!");
          }

          // Check the specific size of each array
          if ((getd(e_np1) != 3) || (gets(e_np1,1) != 3) || (gets(e_np1, 2) !=3)) 
            throw std::runtime_error("e_np1 does not have the right shape!");

          if ((getd(e_n) != 3) || (gets(e_n,1) != 3) || (gets(e_n, 2) !=3)) 
            throw std::runtime_error("e_n does not have the right shape!");

          if ((getd(T_np1) != 1)) 
            throw std::runtime_error("T_np1 does not have the right shape!");

          if ((getd(T_n) != 1)) 
            throw std::runtime_error("T_n does not have the right shape!");

          if ((getd(s_np1) != 3) || (gets(s_np1,1) != 3) || (gets(s_np1, 2) !=3)) 
            throw std::runtime_error("s_np1 does not have the right shape!");

          if ((getd(s_n) != 3) || (gets(s_n,1) != 3) || (gets(s_n, 2) !=3)) 
            throw std::runtime_error("s_n does not have the right shape!");

          if ((getd(h_np1) != 2) || (gets(h_np1,1) != nhist)) 
            throw std::runtime_error("h_np1 does not have the right shape!");

          if ((getd(h_n) != 2) || (gets(h_n,1) != nhist)) 
            throw std::runtime_error("h_n does not have the right shape!");

          if ((getd(A_np1) != 5) || (gets(A_np1,1) != 3) || 
              (gets(A_np1, 2) != 3) || (gets(A_np1,3) != 3) ||
              (gets(A_np1, 4) != 3)) 
            throw std::runtime_error("A_np1 does not have the right shape!");

          if ((getd(u_np1) != 1)) 
            throw std::runtime_error("u_np1 does not have the right shape!");

          if ((getd(u_n) != 1)) 
            throw std::runtime_error("u_n does not have the right shape!");

          if ((getd(p_np1) != 1)) 
            throw std::runtime_error("u_np1 does not have the right shape!");

          if ((getd(p_n) != 1)) 
            throw std::runtime_error("u_n does not have the right shape!");
          
          // Needs to happen before GIL release
          double * e_np1_ptr = arr2ptr<double>(e_np1);
          double * e_n_ptr = arr2ptr<double>(e_n);
          double * T_np1_ptr = arr2ptr<double>(T_np1);
          double * T_n_ptr = arr2ptr<double>(T_n);
          double * s_np1_ptr = arr2ptr<double>(s_np1);
          double * s_n_ptr = arr2ptr<double>(s_n);
          double * h_np1_ptr = arr2ptr<double>(h_np1);
          double * h_n_ptr = arr2ptr<double>(h_n);
          double * A_np1_ptr = arr2ptr<double>(A_np1);
          double * u_np1_ptr = arr2ptr<double>(u_np1);
          double * u_n_ptr = arr2ptr<double>(u_n);
          double * p_np1_ptr = arr2ptr<double>(p_np1);
          double * p_n_ptr = arr2ptr<double>(p_n);

          // Release the GIL and call the batch update
          {
            py::gil_scoped_release release;

            block_evaluate(model, nblock, 
                                 e_np1_ptr, e_n_ptr, T_np1_ptr, T_n_ptr,
                                 t_np1, t_n, s_np1_ptr, s_n_ptr,
                                 h_np1_ptr, h_n_ptr, A_np1_ptr,
                                 u_np1_ptr, u_n_ptr, p_np1_ptr, p_n_ptr);
          }

        }, "Block evaluate a bunch of models in tensor notation");

  m.def("t2m", 
        [](py::array_t<double, py::array::c_style> T) -> py::array_t<double>
        {
          py::buffer_info info = T.request();
          size_t n = info.shape[0];
          auto f = alloc_mat<double>(n, 6);
          t2m(arr2ptr<double>(T), arr2ptr<double>(f), n);
          return f;
        }, "Block convert symmetric tensors to Mandel");
  m.def("m2t",
        [](py::array_t<double, py::array::c_style> T) -> py::array_t<double>
        {
          size_t n = T.request().shape[0];
          auto f = alloc_3d<double>(n, 3, 3);
          m2t(arr2ptr<double>(T), arr2ptr<double>(f), n);
          return f;
        }, "Block convert Mandel vectors to symmetric tensors");
  m.def("m42t4",
        [](py::array_t<double, py::array::c_style> T) -> py::array_t<double>
        {
          size_t n = T.request().shape[0];
          size_t i = 3;
          size_t dim = 5;
          size_t dz = sizeof(double);

          auto f = py::array(py::buffer_info(
                  nullptr, sizeof(double), py::format_descriptor<double>::value,
                  dim, {n,i,i,i,i},
                  {i*i*i*i*dz, i*i*i*dz, i*i*dz, i*dz, dz}));

          m42t4(arr2ptr<double>(T), arr2ptr<double>(f), n);

          return f;
        }, "Block convert Mandel rank 4 tensors to symmetric full tensors");


} // PYBIND11_MODULE

} // namespace neml
