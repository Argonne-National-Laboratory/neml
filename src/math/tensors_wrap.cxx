#include "../pyhelp.h"

#include "tensors.h"

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

              ss << "SymSym(array([";
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
      .def("dot", [](const RankFour & me, const SymSym & other) -> RankFour
           {
            return me.dot(other);
           })
      .def("dot", [](const RankFour & me, const SymSkew & other) -> RankFour
           {
            return me.dot(other);
           })
      .def("dot", [](const RankFour & me, const SkewSym & other) -> RankFour
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
      .def(py::self * SymSym())
      .def(py::self * SymSkew())
      .def(py::self * SkewSym())

      .def(py::self * RankTwo())
      .def(py::self * Symmetric())
      .def(py::self * Skew())
      ;

  py::class_<SymSym, Tensor, std::shared_ptr<SymSym>>(m, "SymSym")
      // Start standard
      .def(py::init<const std::vector<std::vector<double>>>(), py::arg("data"))

      .def("__repr__",
           [](SymSym & me) -> std::string
           {
              std::ostringstream ss;

              ss << "SymSym(array([";
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
           [](SymSym & me) -> std::string
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

      .def("opposite", &SymSym::opposite)
      .def("__neg__", &SymSym::opposite)

      .def(double() * py::self)
      .def(py::self * double())

      .def(py::self / double())

      .def(py::self += py::self)
      .def(py::self + py::self)

      .def(py::self -= py::self)
      .def(py::self - py::self)

      // End standard
      .def("to_full", &SymSym::to_full)

      .def("dot", [](const SymSym & me, const SymSym & other) -> SymSym
           {
            return me.dot(other);
           })
      .def("dot", [](const SymSym & me, const RankFour & other) -> RankFour
           {
            return me.dot(other);
           })
      .def("dot", [](const SymSym & me, const SymSkew & other) -> RankFour
           {
            return me.dot(other);
           })
      .def("dot", [](const SymSym & me, const SymSkew & other) -> RankFour
           {
            return me.dot(other);
           })
        
      .def("dot", [](const SymSym & me, const Symmetric & other) -> Symmetric
           {
            return me.dot(other);
           })
      .def("dot", [](const SymSym & me, const RankTwo & other) -> RankTwo
           {
            return me.dot(other);
           })
      .def("dot", [](const SymSym & me, const Skew & other) -> RankTwo
           {
            return me.dot(other);
           })

      .def(py::self * py::self)
      .def(py::self * RankFour())
      .def(py::self * SymSkew())
      .def(py::self * SkewSym())
      
      .def(py::self * Symmetric())
      .def(py::self * RankTwo())
      .def(py::self * Skew())
      ;

  py::class_<SymSkew, Tensor, std::shared_ptr<SymSkew>>(m, "SymSkew")
      // Start standard
      .def(py::init<const std::vector<std::vector<double>>>(), py::arg("data"))

      .def("__repr__",
           [](SymSkew & me) -> std::string
           {
              std::ostringstream ss;

              ss << "SymSkew(array([";
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
           [](SymSkew & me) -> std::string
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

      .def("opposite", &SymSkew::opposite)
      .def("__neg__", &SymSkew::opposite)

      .def(double() * py::self)
      .def(py::self * double())

      .def(py::self / double())

      .def(py::self += py::self)
      .def(py::self + py::self)

      .def(py::self -= py::self)
      .def(py::self - py::self)

      // End standard
      .def("to_full", &SymSkew::to_full)

      .def("dot", [](const SymSkew & me, const RankFour & other) -> RankFour
           {
            return me.dot(other);
           })
      .def("dot", [](const SymSkew & me, const SymSkew & other) -> RankFour
           {
            return me.dot(other);
           })
      .def("dot", [](const SymSkew & me, const SymSym & other) -> RankFour
           {
            return me.dot(other);
           })
      .def("dot", [](const SymSkew & me, const SkewSym & other) -> RankFour
           {
            return me.dot(other);
           })

      .def("dot", [](const SymSkew & me, const RankTwo & other) -> RankTwo
           {
            return me.dot(other);
           })
      .def("dot", [](const SymSkew & me, const Symmetric & other) -> RankTwo
           {
            return me.dot(other);
           })
      .def("dot", [](const SymSkew & me, const Skew & other) -> RankTwo
           {
            return me.dot(other);
           })

      .def(py::self * py::self)
      .def(py::self * SymSym())
      .def(py::self * RankFour())
      .def(py::self * SkewSym())

      .def(py::self * RankTwo())
      .def(py::self * Symmetric())
      .def(py::self * Skew())
      ;

  py::class_<SkewSym, Tensor, std::shared_ptr<SkewSym>>(m, "SkewSym")
      // Start standard
      .def(py::init<const std::vector<std::vector<double>>>(), py::arg("data"))

      .def("__repr__",
           [](SkewSym & me) -> std::string
           {
              std::ostringstream ss;

              ss << "SkewSym(array([";
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
           [](SkewSym & me) -> std::string
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

      .def("opposite", &SkewSym::opposite)
      .def("__neg__", &SkewSym::opposite)

      .def(double() * py::self)
      .def(py::self * double())

      .def(py::self / double())

      .def(py::self += py::self)
      .def(py::self + py::self)

      .def(py::self -= py::self)
      .def(py::self - py::self)

      // End standard
      .def("to_full", &SkewSym::to_full)

      .def("dot", [](const SkewSym & me, const RankFour & other) -> RankFour
           {
            return me.dot(other);
           })
      .def("dot", [](const SkewSym & me, const SkewSym & other) -> RankFour
           {
            return me.dot(other);
           })
      .def("dot", [](const SkewSym & me, const SymSym & other) -> RankFour
           {
            return me.dot(other);
           })
      .def("dot", [](const SkewSym & me, const SymSkew & other) -> RankFour
           {
            return me.dot(other);
           })

      .def("dot", [](const SkewSym & me, const RankTwo & other) -> RankTwo
           {
            return me.dot(other);
           })
      .def("dot", [](const SkewSym & me, const Symmetric & other) -> RankTwo
           {
            return me.dot(other);
           })
      .def("dot", [](const SkewSym & me, const Skew & other) -> RankTwo
           {
            return me.dot(other);
           })

      .def(py::self * py::self)
      .def(py::self * SymSym())
      .def(py::self * RankFour())
      .def(py::self * SymSkew())

      .def(py::self * RankTwo())
      .def(py::self * Symmetric())
      .def(py::self * Skew())
      ;

  m.def("douter", [](const Symmetric & a, const Symmetric & b) -> SymSym
        {
          return douter(a, b);
        });
  m.def("douter", [](const Skew & a, const Symmetric & b) -> SkewSym
        {
          return douter(a, b);
        });
  m.def("SymSymSkew_SkewSymSym", &SymSymSkew_SkewSymSym);
  m.def("SymSkewSym_SkewSymSym", &SymSkewSym_SkewSymSym);
  m.def("SpecialSymSymSym", &SpecialSymSymSym);

} // PYBIND11_MODULE(tensors, m)

} // namespace neml
