### Macro helper for python bindings ###
macro(pybind sname)
      set(wrapper "${sname}_wrap.cxx")
      pybind11_add_module(${sname} SHARED ${wrapper})
      target_link_libraries(${sname} PRIVATE neml)
endmacro(pybind)
