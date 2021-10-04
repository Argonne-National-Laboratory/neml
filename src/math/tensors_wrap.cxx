#include "pyhelp.h"

#include "math/tensors.h"

#include <iostream>
#include <string>
#include <stdexcept>

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml {

PYBIND11_MODULE(tensors, m) {
  m.doc() = "Standard tensor properties";

  py::class_<Tensor, std::shared_ptr<Tensor>>(m, "Tensor")
      .def_property_readonly("istore", &Tensor::istore)
      .def_property_readonly("data",
           [](Tensor & me) -> py::array_t<double>
           {
            auto v = alloc_vec<double>(me.n());
            std::copy(me.data(), me.data()+me.n(), arr2ptr<double>(v));
            return v;
           }, "raw data view")

      .def(py::self *= double())
      .def(py::self /= double())

      .def(py::self == py::self)
      .def(py::self != py::self)
      ;

  py::class_<Vector, Tensor, std::shared_ptr<Vector>>(m, "Vector")
      // Start standard
      .def(py::init<const std::vector<double>>(), py::arg("data"))

      .def("__repr__",
           [](Vector & me) -> std::string
           {
              std::ostringstream ss;
              ss << "Vector(array([" << me.data()[0] << ", "
                << me.data()[1] << ", " << me.data()[2] << ", " << "]))";
              return ss.str();
           }, "python __repr__")

      .def("__str__",
           [](Vector & me) -> std::string
           {
              std::ostringstream ss;
              ss << "[" << me.data()[0] << ", "
                << me.data()[1] << ", " << me.data()[2] << ", " << "]";
              return ss.str();
           }, "python __str__")


      .def("opposite", &Vector::opposite)
      .def("__neg__", &Vector::opposite)

      .def(double() * py::self)
      .def(py::self * double())

      .def(py::self / double())

      .def(py::self += py::self)
      .def(py::self + py::self)

      .def(py::self -= py::self)
      .def(py::self - py::self)

      .def("__getitem__", [](const Vector & v, size_t i) {
           if (i >= v.n()) throw py::index_error();
           return v(i);
      })

      .def("__setitem__", [](Vector & v, size_t i, double val) {
           if (i >= v.n()) throw py::index_error();
           v(i) = val;
      })

      // End standard

      .def("dot", &Vector::dot)
      .def("outer", &Vector::outer)
      .def("norm", &Vector::norm)
      .def("cross", &Vector::cross)
      .def("normalize", &Vector::normalize)
      ;

  m.def("outer", &outer);

  py::class_<RankTwo, Tensor, std::shared_ptr<RankTwo>>(m, "RankTwo")
      // Start standard
      .def(py::init<const std::vector<std::vector<double>>>(), py::arg("data"))

      .def("__repr__",
           [](RankTwo & me) -> std::string
           {
              std::ostringstream ss;
              
              ss << "RankTwo(array([";
              for (size_t i=0; i<3; i++) {
                ss << "[";

                for (size_t j=0; j<3; j++) {
                  ss << me.data()[i*3+j] << " ";
                }
                ss << "]" << std::endl;
              }
              ss << "]))";

              return ss.str();
           }, "python __repr__")

      .def("__str__",
           [](RankTwo & me) -> std::string
           {
              std::ostringstream ss;
              
              ss << "[";
              for (size_t i=0; i<3; i++) {
                ss << "[";

                for (size_t j=0; j<3; j++) {
                  ss << me.data()[i*3+j] << " ";
                }
                ss << "]" << std::endl;
              }
              ss << "]";

              return ss.str();
           }, "python __str__")

      .def("opposite", &RankTwo::opposite)
      .def("__neg__", &RankTwo::opposite)

      .def(double() * py::self)
      .def(py::self * double())

      .def(py::self / double())

      .def(py::self += py::self)
      .def(py::self + py::self)
      .def(py::self += Symmetric())
      .def(py::self + Symmetric())
      .def(Symmetric() + py::self)
      .def(py::self += Skew())
      .def(py::self + Skew())
      .def(Skew() + py::self)

      .def(py::self -= py::self)
      .def(py::self - py::self)
      .def(py::self -= Symmetric())
      .def(py::self - Symmetric())
      .def(Symmetric() - py::self)
      .def(py::self -= Skew())
      .def(py::self - Skew())
      .def(Skew() - py::self)

      .def_static("id", &Symmetric::id)

      .def("__getitem__", [](const RankTwo & M, std::tuple<size_t,size_t> ind) {
           size_t i = std::get<0>(ind);
           size_t j = std::get<1>(ind);
           if ((i >= 3) || (j > 3)) throw py::index_error();
           return M(i,j);
      })

      .def("__setitem__", [](RankTwo & M, std::tuple<size_t,size_t> ind, double val) {
           size_t i = std::get<0>(ind);
           size_t j = std::get<1>(ind);
           if ((i >= 3) || (j > 3)) throw py::index_error();
           M(i,j) = val;
      })

      // End standard
      .def("dot", [](const RankTwo & me, const Vector & other) -> Vector
           {
            return me.dot(other);
           })
      .def("dot", [](const RankTwo & me, const RankTwo & other) -> RankTwo
           {
            return me.dot(other);
           })
      .def("dot", [](const RankTwo & me, const Symmetric & other) -> RankTwo
           {
            return me.dot(other);
           })
      .def("dot", [](const RankTwo & me, const Skew & other) -> RankTwo
           {
            return me.dot(other);
           })

      .def("contract", (double (RankTwo::*)(const RankTwo &) const) 
           &RankTwo::contract)
      .def("contract", (double (RankTwo::*)(const Symmetric &) const) 
           &RankTwo::contract)
      .def("contract", (double (RankTwo::*)(const Skew &) const) 
           &RankTwo::contract)

      .def(py::self * py::self)
      .def(py::self * Vector())
      .def(Vector() * py::self)
      .def(py::self * Symmetric())
      .def(Symmetric() * py::self)
      .def(py::self * Skew())
      .def(Skew() * py::self)

      .def("inverse", &RankTwo::inverse)
      .def("transpose", &RankTwo::transpose)
      .def("norm", &RankTwo::norm)
      ;

  py::class_<Symmetric, Tensor, std::shared_ptr<Symmetric>>(m, "Symmetric")
      // Start standard
      .def(py::init<const RankTwo &>(), py::arg("Full"))
      .def(py::init<const std::vector<std::vector<double>>>(), py::arg("data"))

      .def("__repr__",
           [](Symmetric & me) -> std::string
           {
              std::ostringstream ss;

              RankTwo you(me);
              
              ss << "RankTwo(array([";
              for (size_t i=0; i<3; i++) {
                ss << "[";

                for (size_t j=0; j<3; j++) {
                  ss << you.data()[i*3+j] << " ";
                }
                ss << "]" << std::endl;
              }
              ss << "]))";

              return ss.str();
           }, "python __repr__")

      .def("__str__",
           [](Symmetric & me) -> std::string
           {
              std::ostringstream ss;

              RankTwo you(me);
              
              ss << "[";
              for (size_t i=0; i<3; i++) {
                ss << "[";

                for (size_t j=0; j<3; j++) {
                  ss << you.data()[i*3+j] << " ";
                }
                ss << "]" << std::endl;
              }
              ss << "]";

              return ss.str();
           }, "python __str__")

      .def("opposite", &Symmetric::opposite)
      .def("__neg__", &Symmetric::opposite)

      .def(double() * py::self)
      .def(py::self * double())

      .def(py::self / double())

      .def(py::self += py::self)
      .def(py::self + py::self)

      .def(py::self -= py::self)
      .def(py::self - py::self)

      // End standard
      
      .def_static("id", &Symmetric::id)
      .def("trace", &Symmetric::trace)
      .def("dev", &Symmetric::dev)

      .def("dot", [](const Symmetric & me, const Vector & other) -> Vector
           {
            return me.dot(other);
           })
      .def("dot", [](const Symmetric & me, const Symmetric & other) -> Symmetric
           {
            return me.dot(other);
           })
      .def("dot", [](const Symmetric & me, const RankTwo & other) -> RankTwo
           {
            return me.dot(other);
           })
      .def("dot", [](const Symmetric & me, const Skew & other) -> RankTwo
           {
            return me.dot(other);
           })

      .def("contract", (double (Symmetric::*)(const RankTwo &) const) 
           &Symmetric::contract)
      .def("contract", (double (Symmetric::*)(const Symmetric &) const) 
           &Symmetric::contract)
      .def("contract", (double (Symmetric::*)(const Skew &) const) 
           &Symmetric::contract)

      .def(py::self * py::self)
      .def(py::self * Vector())
      .def(Vector() * py::self)

      .def(py::self * Skew())

      .def("inverse", &Symmetric::inverse)
      .def("transpose", &Symmetric::transpose)
      .def("norm", &Symmetric::norm)
      .def("to_full", &Symmetric::to_full)
      ;

  py::class_<Skew, Tensor, std::shared_ptr<Skew>>(m, "Skew")
      // Start standard
      .def(py::init<const RankTwo &>(), py::arg("Full"))
      .def(py::init<const std::vector<std::vector<double>>>(), py::arg("data"))

      .def("__repr__",
           [](Skew & me) -> std::string
           {
              std::ostringstream ss;

              RankTwo you(me);
              
              ss << "RankTwo(array([";
              for (size_t i=0; i<3; i++) {
                ss << "[";

                for (size_t j=0; j<3; j++) {
                  ss << you.data()[i*3+j] << " ";
                }
                ss << "]" << std::endl;
              }
              ss << "]))";

              return ss.str();
           }, "python __repr__")

      .def("__str__",
           [](Skew & me) -> std::string
           {
              std::ostringstream ss;

              RankTwo you(me);
              
              ss << "[";
              for (size_t i=0; i<3; i++) {
                ss << "[";

                for (size_t j=0; j<3; j++) {
                  ss << you.data()[i*3+j] << " ";
                }
                ss << "]" << std::endl;
              }
              ss << "]";

              return ss.str();
           }, "python __str__")

      .def("opposite", &Skew::opposite)
      .def("__neg__", &Skew::opposite)

      .def(double() * py::self)
      .def(py::self * double())

      .def(py::self / double())

      .def(py::self += py::self)
      .def(py::self + py::self)

      .def(py::self -= py::self)
      .def(py::self - py::self)

      // End standard
      .def("dot", [](const Skew & me, const Vector & other) -> Vector
           {
            return me.dot(other);
           })
      .def("dot", [](const Skew & me, const Skew & other) -> Skew
           {
            return me.dot(other);
           })
      .def("dot", [](const Skew & me, const Symmetric & other) -> RankTwo
           {
            return me.dot(other);
           })
      .def("dot", [](const Skew & me, const RankTwo & other) -> RankTwo
           {
            return me.dot(other);
           })

      .def("contract", (double (Skew::*)(const RankTwo &) const) 
           &Skew::contract)
      .def("contract", (double (Skew::*)(const Symmetric &) const) 
           &Skew::contract)
      .def("contract", (double (Skew::*)(const Skew &) const) 
           &Skew::contract)

      .def(py::self * py::self)
      .def(py::self * Vector())
      .def(Vector() * py::self)

      .def(py::self * Symmetric())

      .def("transpose", &Skew::transpose)
      ;

  py::class_<RankFour, Tensor, std::shared_ptr<RankFour>>(m, "RankFour")
      // Start standard
      .def(py::init<const std::vector<std::vector<std::vector<std::vector<double>>>>>(), py::arg("data"))

      .def("__repr__",
           [](RankFour & me) -> std::string
           {
              std::ostringstream ss;

              ss << "RankFour(array([";
              for (size_t i=0; i<9; i++) {
                ss << "[";

                for (size_t j=0; j<9; j++) {
                  ss << me.data()[i*9+j] << " ";
                }
                ss << "]" << std::endl;
              }
              ss << "]))";

              return ss.str();
           }, "python __repr__")

      .def("__str__",
           [](RankFour & me) -> std::string
           {
              std::ostringstream ss;

              ss << "[";
              for (size_t i=0; i<9; i++) {
                ss << "[";

                for (size_t j=0; j<9; j++) {
                  ss << me.data()[i*9+j] << " ";
                }
                ss << "]" << std::endl;
              }
              ss << "]";

              return ss.str();
           }, "python __str__")

      .def("opposite", &RankFour::opposite)
      .def("__neg__", &RankFour::opposite)

      .def(double() * py::self)
      .def(py::self * double())

      .def(py::self / double())

      .def(py::self += py::self)
      .def(py::self + py::self)

      .def(py::self -= py::self)
      .def(py::self - py::self)

      // End standard
      .def("to_sym", &RankFour::to_sym)
      .def("to_symskew", &RankFour::to_symskew)
      .def("to_skewsym", &RankFour::to_skewsym)

      .def("__getitem__", [](const RankFour & M, std::tuple<size_t,size_t,size_t,size_t> ind) {
           size_t i = std::get<0>(ind);
           size_t j = std::get<1>(ind);
           size_t k = std::get<2>(ind);
           size_t l = std::get<3>(ind);
           if ((i >= 3) || (j > 3) || (k>3) || (l>3)) throw py::index_error();
           return M(i,j,k,l);
      })

      .def("__setitem__", [](RankFour & M, std::tuple<size_t,size_t,size_t,size_t> ind, double val) {
           size_t i = std::get<0>(ind);
           size_t j = std::get<1>(ind);
           size_t k = std::get<2>(ind);
           size_t l = std::get<3>(ind);
           if ((i >= 3) || (j > 3) || (k>3) || (l>3)) throw py::index_error();
           M(i,j,k,l) = val;
      })
            
      .def("dot", [](const RankFour & me, const RankFour & other) -> RankFour
           {
            return me.dot(other);
           })
      .def("dot", [](const RankFour & me, const SymSymR4 & other) -> RankFour
           {
            return me.dot(other);
           })
      .def("dot", [](const RankFour & me, const SymSkewR4 & other) -> RankFour
           {
            return me.dot(other);
           })
      .def("dot", [](const RankFour & me, const SkewSymR4 & other) -> RankFour
           {
            return me.dot(other);
          })

      .def("dot", [](const RankFour & me, const RankTwo & other) -> RankTwo
           {
            return me.dot(other);
           })
      .def("dot", [](const RankFour & me, const Symmetric & other) -> RankTwo
           {
            return me.dot(other);
           })
      .def("dot", [](const RankFour & me, const Skew & other) -> RankTwo
           {
            return me.dot(other);
           })

      .def(py::self * py::self)
      .def(py::self * SymSymR4())
      .def(py::self * SymSkewR4())
      .def(py::self * SkewSymR4())

      .def(py::self * RankTwo())
      .def(py::self * Symmetric())
      .def(py::self * Skew())
      ;

  py::class_<SymSymR4, Tensor, std::shared_ptr<SymSymR4>>(m, "SymSymR4")
      // Start standard
      .def(py::init<const std::vector<std::vector<double>>>(), py::arg("data"))

      .def("__repr__",
           [](SymSymR4 & me) -> std::string
           {
              std::ostringstream ss;

              ss << "SymSymR4(array([";
              for (size_t i=0; i<6; i++) {
                ss << "[";

                for (size_t j=0; j<6; j++) {
                  ss << me.data()[i*6+j] << " ";
                }
                ss << "]" << std::endl;
              }
              ss << "]))";

              return ss.str();
           }, "python __repr__")

      .def("__str__",
           [](SymSymR4 & me) -> std::string
           {
              std::ostringstream ss;

              ss << "[";
              for (size_t i=0; i<6; i++) {
                ss << "[";

                for (size_t j=0; j<6; j++) {
                  ss << me.data()[i*6+j] << " ";
                }
                ss << "]" << std::endl;
              }
              ss << "]";

              return ss.str();
           }, "python __str__")

      .def("opposite", &SymSymR4::opposite)
      .def("__neg__", &SymSymR4::opposite)

      .def(double() * py::self)
      .def(py::self * double())

      .def(py::self / double())

      .def(py::self += py::self)
      .def(py::self + py::self)

      .def(py::self -= py::self)
      .def(py::self - py::self)

      // End standard
      .def("to_full", &SymSymR4::to_full)

      .def("dot", [](const SymSymR4 & me, const SymSymR4 & other) -> SymSymR4
           {
            return me.dot(other);
           })
      .def("dot", [](const SymSymR4 & me, const RankFour & other) -> RankFour
           {
            return me.dot(other);
           })
      .def("dot", [](const SymSymR4 & me, const SymSkewR4 & other) -> RankFour
           {
            return me.dot(other);
           })
      .def("dot", [](const SymSymR4 & me, const SymSkewR4 & other) -> RankFour
           {
            return me.dot(other);
           })
        
      .def("dot", [](const SymSymR4 & me, const Symmetric & other) -> Symmetric
           {
            return me.dot(other);
           })
      .def("dot", [](const SymSymR4 & me, const RankTwo & other) -> RankTwo
           {
            return me.dot(other);
           })
      .def("dot", [](const SymSymR4 & me, const Skew & other) -> RankTwo
           {
            return me.dot(other);
           })

      .def(py::self * py::self)
      .def(py::self * RankFour())
      .def(py::self * SymSkewR4())
      .def(py::self * SkewSymR4())
      
      .def(py::self * Symmetric())
      .def(py::self * RankTwo())
      .def(py::self * Skew())

      .def_static("id", &SymSymR4::id)
      .def_static("id_dev", &SymSymR4::id_dev)
      ;

  py::class_<SymSkewR4, Tensor, std::shared_ptr<SymSkewR4>>(m, "SymSkewR4")
      // Start standard
      .def(py::init<const std::vector<std::vector<double>>>(), py::arg("data"))

      .def("__repr__",
           [](SymSkewR4 & me) -> std::string
           {
              std::ostringstream ss;

              ss << "SymSkewR4(array([";
              for (size_t i=0; i<6; i++) {
                ss << "[";

                for (size_t j=0; j<3; j++) {
                  ss << me.data()[i*3+j] << " ";
                }
                ss << "]" << std::endl;
              }
              ss << "]))";

              return ss.str();
           }, "python __repr__")

      .def("__str__",
           [](SymSkewR4 & me) -> std::string
           {
              std::ostringstream ss;

              ss << "[";
              for (size_t i=0; i<6; i++) {
                ss << "[";

                for (size_t j=0; j<3; j++) {
                  ss << me.data()[i*3+j] << " ";
                }
                ss << "]" << std::endl;
              }
              ss << "]";

              return ss.str();
           }, "python __str__")

      .def("opposite", &SymSkewR4::opposite)
      .def("__neg__", &SymSkewR4::opposite)

      .def(double() * py::self)
      .def(py::self * double())

      .def(py::self / double())

      .def(py::self += py::self)
      .def(py::self + py::self)

      .def(py::self -= py::self)
      .def(py::self - py::self)

      // End standard
      .def("to_full", &SymSkewR4::to_full)

      .def("dot", [](const SymSkewR4 & me, const RankFour & other) -> RankFour
           {
            return me.dot(other);
           })
      .def("dot", [](const SymSkewR4 & me, const SymSkewR4 & other) -> RankFour
           {
            return me.dot(other);
           })
      .def("dot", [](const SymSkewR4 & me, const SymSymR4 & other) -> RankFour
           {
            return me.dot(other);
           })
      .def("dot", [](const SymSkewR4 & me, const SkewSymR4 & other) -> RankFour
           {
            return me.dot(other);
           })

      .def("dot", [](const SymSkewR4 & me, const RankTwo & other) -> RankTwo
           {
            return me.dot(other);
           })
      .def("dot", [](const SymSkewR4 & me, const Symmetric & other) -> RankTwo
           {
            return me.dot(other);
           })
      .def("dot", [](const SymSkewR4 & me, const Skew & other) -> RankTwo
           {
            return me.dot(other);
           })

      .def(py::self * py::self)
      .def(py::self * SymSymR4())
      .def(py::self * RankFour())
      .def(py::self * SkewSymR4())

      .def(py::self * RankTwo())
      .def(py::self * Symmetric())
      .def(py::self * Skew())
      ;

  py::class_<SkewSymR4, Tensor, std::shared_ptr<SkewSymR4>>(m, "SkewSymR4")
      // Start standard
      .def(py::init<const std::vector<std::vector<double>>>(), py::arg("data"))

      .def("__repr__",
           [](SkewSymR4 & me) -> std::string
           {
              std::ostringstream ss;

              ss << "SkewSymR4(array([";
              for (size_t i=0; i<3; i++) {
                ss << "[";

                for (size_t j=0; j<6; j++) {
                  ss << me.data()[i*6+j] << " ";
                }
                ss << "]" << std::endl;
              }
              ss << "]))";

              return ss.str();
           }, "python __repr__")

      .def("__str__",
           [](SkewSymR4 & me) -> std::string
           {
              std::ostringstream ss;

              ss << "[";
              for (size_t i=0; i<3; i++) {
                ss << "[";

                for (size_t j=0; j<6; j++) {
                  ss << me.data()[i*6+j] << " ";
                }
                ss << "]" << std::endl;
              }
              ss << "]";

              return ss.str();
           }, "python __str__")

      .def("opposite", &SkewSymR4::opposite)
      .def("__neg__", &SkewSymR4::opposite)

      .def(double() * py::self)
      .def(py::self * double())

      .def(py::self / double())

      .def(py::self += py::self)
      .def(py::self + py::self)

      .def(py::self -= py::self)
      .def(py::self - py::self)

      // End standard
      .def("to_full", &SkewSymR4::to_full)

      .def("dot", [](const SkewSymR4 & me, const RankFour & other) -> RankFour
           {
            return me.dot(other);
           })
      .def("dot", [](const SkewSymR4 & me, const SkewSymR4 & other) -> RankFour
           {
            return me.dot(other);
           })
      .def("dot", [](const SkewSymR4 & me, const SymSymR4 & other) -> RankFour
           {
            return me.dot(other);
           })
      .def("dot", [](const SkewSymR4 & me, const SymSkewR4 & other) -> RankFour
           {
            return me.dot(other);
           })

      .def("dot", [](const SkewSymR4 & me, const RankTwo & other) -> RankTwo
           {
            return me.dot(other);
           })
      .def("dot", [](const SkewSymR4 & me, const Symmetric & other) -> RankTwo
           {
            return me.dot(other);
           })
      .def("dot", [](const SkewSymR4 & me, const Skew & other) -> RankTwo
           {
            return me.dot(other);
           })

      .def(py::self * py::self)
      .def(py::self * SymSymR4())
      .def(py::self * RankFour())
      .def(py::self * SymSkewR4())

      .def(py::self * RankTwo())
      .def(py::self * Symmetric())
      .def(py::self * Skew())
      ;

  m.def("douter", [](const Symmetric & a, const Symmetric & b) -> SymSymR4
        {
          return douter(a, b);
        });
  m.def("douter", [](const Skew & a, const Symmetric & b) -> SkewSymR4
        {
          return douter(a, b);
        });
  m.def("SymSymR4Skew_SkewSymR4SymR4", &SymSymR4Skew_SkewSymR4SymR4);
  m.def("SymSkewR4Sym_SkewSymR4SymR4", &SymSkewR4Sym_SkewSymR4SymR4);
  m.def("SpecialSymSymR4Sym", &SpecialSymSymR4Sym);

  py::class_<SymSymSymR6, Tensor, std::shared_ptr<SymSymSymR6>>(m, "SymSymSymR6")
      // Start standard
      .def(py::init<const std::vector<std::vector<std::vector<double>>>>(), py::arg("data"))

      .def("__repr__",
           [](SymSymSymR6 & me) -> std::string
           {
              std::ostringstream ss;

              ss << "SymSymSymR6(array([";
              for (size_t i=0; i<6; i++) {
                ss << "[";

                for (size_t j=0; j<6; j++) {
                ss << "[";
                  for (size_t k=0; k<6; k++) {
                    ss << me.data()[i*36+j*6+k] << " ";
                  }
                  ss << "]" << std::endl;
                }
                ss << "]" << std::endl;
              }
              ss << "]))";

              return ss.str();
           }, "python __repr__")

      .def("__str__",
           [](SymSymSymR6 & me) -> std::string
           {
              std::ostringstream ss;

              ss << "[";
              for (size_t i=0; i<6; i++) {
                ss << "[";

                for (size_t j=0; j<6; j++) {
                  ss << "[";
                  for (size_t k=0; k<6; k++) {
                    ss << me.data()[i*36+j*6+k] << " ";
                  }
                  ss << "]" << std::endl;
                }
                ss << "]" << std::endl;
              }
              ss << "]";

              return ss.str();
           }, "python __str__")

      .def("opposite", &SymSymSymR6::opposite)
      .def("__neg__", &SymSymSymR6::opposite)

      .def(double() * py::self)
      .def(py::self * double())

      .def(py::self / double())

      .def(py::self += py::self)
      .def(py::self + py::self)

      .def(py::self -= py::self)
      .def(py::self - py::self)

      .def("dot_i", &SymSymSymR6::dot_i)
      .def("dot_j", &SymSymSymR6::dot_j)
      .def("dot_k", &SymSymSymR6::dot_k)
      ;

      m.def("outer_product_k", &outer_product_k);

} // PYBIND11_MODULE(tensors, m)

} // namespace neml
