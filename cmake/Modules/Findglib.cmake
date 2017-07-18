# Variables to define:
#  Glib_FOUND - system has Glib
#  Glib_INCLUDE_DIRS - the Glib include directories
#  Glib_LIBRARIES - link these to use Glib

# Hunt down pkg-config info
find_package(PkgConfig)
pkg_check_modules(Glib_PKGCONF glib-2.0 REQUIRED)

# Main include dir
find_path(Glib_INCLUDE_DIR
  NAMES glib.h
  PATHS ${Glib_PKGCONF_INCLUDE_DIRS}
  PATH_SUFFIXES glib-2.0
)

# Glib-related libraries also use a separate config header, which is in lib dir
find_path(GlibConfig_INCLUDE_DIR
  NAMES glibconfig.h
  PATHS ${Glib_PKGCONF_INCLUDE_DIRS} /usr
  PATH_SUFFIXES lib/glib-2.0/include
)

# Finally the library itself
find_library(Glib_LIBRARY
  NAMES glib-2.0
  PATHS ${Glib_PKGCONF_LIBRARY_DIRS}
)

set(Glib_FOUND TRUE)
set(Glib_INCLUDE_DIRS ${Glib_INCLUDE_DIR} ${GlibConfig_INCLUDE_DIR})
set(Glib_LIBRARIES ${Glib_LIBRARY})

mark_as_advanced(GlibConfig_INCLUDE_DIR Glib_INCLUDE_DIR Glib_LIBRARY)
