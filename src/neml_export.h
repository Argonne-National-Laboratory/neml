
#ifndef NEML_EXPORT_H
#define NEML_EXPORT_H

#ifdef NEML_STATIC_DEFINE
#  define NEML_EXPORT
#  define NEML_NO_EXPORT
#else
#  ifndef NEML_EXPORT
#    ifdef neml_EXPORTS
        /* We are building this library */
#      define NEML_EXPORT __declspec(dllexport)
#    else
        /* We are using this library */
#      define NEML_EXPORT __declspec(dllimport)
#    endif
#  endif

#  ifndef NEML_NO_EXPORT
#    define NEML_NO_EXPORT 
#  endif
#endif

#ifndef NEML_DEPRECATED
#  define NEML_DEPRECATED __declspec(deprecated)
#endif

#ifndef NEML_DEPRECATED_EXPORT
#  define NEML_DEPRECATED_EXPORT NEML_EXPORT NEML_DEPRECATED
#endif

#ifndef NEML_DEPRECATED_NO_EXPORT
#  define NEML_DEPRECATED_NO_EXPORT NEML_NO_EXPORT NEML_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef NEML_NO_DEPRECATED
#    define NEML_NO_DEPRECATED
#  endif
#endif

#endif /* NEML_EXPORT_H */
