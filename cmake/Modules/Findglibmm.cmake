# Variables to define:
#  glibmm_FOUND - system has Glib
#  glibmm_INCLUDE_DIRS - the Glib include directories
#  glibmm_LIBRARIES - link these to use Glib

IF(glibmm_FIND_REQUIRED)
    FIND_PACKAGE(glib REQUIRED)
    FIND_PACKAGE(sigc++ REQUIRED)
ELSE(glibmm_FIND_REQUIRED)
    FIND_PACKAGE(glib)
    FIND_PACKAGE(sigc++)
ENDIF(glibmm_FIND_REQUIRED)

# Hunt down pkg-config info
find_package(PkgConfig)
pkg_check_modules(glibmm_PKGCONF glibmm-2.4 REQUIRED)

# Main include dir
find_path(glibmm_INCLUDE_DIR
  NAMES glibmm/main.h
  PATHS ${glibmm_PKGCONF_INCLUDE_DIRS}
  PATH_SUFFIXES glibmm-2.4
)

# Glib-related libraries also use a separate config header, which is in lib dir
find_path(glibmmConfig_INCLUDE_DIR
  NAMES glibmmconfig.h
  PATHS ${glibmm_PKGCONF_INCLUDE_DIRS}
  PATH_SUFFIXES lib/glibmm-2.4/include
)

# Finally the library itself
find_library(glibmm_LIBRARY
  NAMES glibmm glibmm-2.4
  PATHS ${glibmm_PKGCONF_LIBRARY_DIRS}
)

set(glibmm_FOUND TRUE)
set(glibmm_INCLUDE_DIRS ${glibmm_INCLUDE_DIR} ${glibmmConfig_INCLUDE_DIR} ${Glib_INCLUDE_DIRS} ${sigc++_INCLUDE_DIRS})
set(glibmm_LIBRARIES ${glibmm_LIBRARY} ${Glib_LIBRARES} ${sigc++_LIBRARIES})

mark_as_advanced(glibmmConfig_INCLUDE_DIR glibmm_INCLUDE_DIR glibmm_LIBRARY)


