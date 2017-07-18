# Variables to define:
#  sigc++_FOUND - system has Glib
#  sigc++_INCLUDE_DIRS - the Glib include directories
#  sigc++_LIBRARIES - link these to use Glib

# Hunt down pkg-config info
find_package(PkgConfig)
pkg_check_modules(sigc++_PKGCONF sigc++-2.0 REQUIRED)

# Main include dir
find_path(sigc++_INCLUDE_DIR
  NAMES sigc++/sigc++.h
  PATHS ${sigc++_PKGCONF_INCLUDE_DIRS}
  PATH_SUFFIXES sigc++-2.0
)

# Glib-related libraries also use a separate config header, which is in lib dir
find_path(sigc++Config_INCLUDE_DIR
  NAMES sigc++config.h
  PATHS ${sigc++_PKGCONF_INCLUDE_DIRS}
  PATH_SUFFIXES lib/sigc++-2.0/include
)

# Finally the library itself
find_library(sigc++_LIBRARY
  NAMES sigc sigc-2.0
  PATHS ${sigc++_PKGCONF_LIBRARY_DIRS}
)

set(sigc++_FOUND TRUE)
set(sigc++_INCLUDE_DIRS ${sigc++_INCLUDE_DIR} ${sigc++Config_INCLUDE_DIR})
set(sigc++_LIBRARIES ${sigc++_LIBRARY})

mark_as_advanced(sigc++Config_INCLUDE_DIR sigc++_INCLUDE_DIR sigc++_LIBRARY)


