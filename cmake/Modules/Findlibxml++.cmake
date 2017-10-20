# Variables to define:
#  libxml++_FOUND - system has Glib
#  libxml++_INCLUDE_DIRS - the Glib include directories
#  libxml++_LIBRARIES - link these to use Glib

IF(libxml++_FIND_REQUIRED)
    FIND_PACKAGE(glibmm REQUIRED)
    FIND_PACKAGE(LibXml2 REQUIRED)
ELSE(libxml++_FIND_REQUIRED)
    FIND_PACKAGE(glibmm)
    FIND_PACKAGE(LibXml2)
ENDIF(libxml++_FIND_REQUIRED)

# Hunt down pkg-config info
find_package(PkgConfig)

pkg_search_module(libxml++_PKGCONF libxml++-3.0)
if(${libxml++_PKGCONF_FOUND})
      add_definitions(-DLIBXML++V3)
else()
      pkg_search_module(libxml++_PKGCONF REQUIRED libxml++-2.9 libxml++-2.8 libxml++-2.7 libxml++-2.6)
endif()


# Main include dir
find_path(libxml++_INCLUDE_DIR
  NAMES libxml++/libxml++.h
  PATHS ${libxml++_PKGCONF_INCLUDE_DIRS}
)

# Glib-related libraries also use a separate config header, which is in lib dir
find_path(libxml++Config_INCLUDE_DIR
  NAMES libxml++config.h
  PATHS ${libxml++_PKGCONF_INCLUDE_DIRS}
)

# Finally the library itself
find_library(libxml++_LIBRARY
  NAMES xml++ xml++-3.0 xml++-2.9 xml++-2.8 xml++-2.7 xml++-2.6
  PATHS ${libxml++_PKGCONF_LIBRARY_DIRS}
)

set(libxml++_FOUND TRUE)
set(libxml++_INCLUDE_DIRS ${libxml++_INCLUDE_DIR} ${libxml++Config_INCLUDE_DIR} ${glibmm_INCLUDE_DIRS} ${LibXml2_INLUDE_DIRS})
set(libxml++_LIBRARIES ${libxml++_LIBRARY} ${glibmm_LIBRARIES} ${LibXml2_LIBRARIES})

mark_as_advanced(libxml++Config_INCLUDE_DIR libxml++_INCLUDE_DIR libxml++_LIBRARY)



