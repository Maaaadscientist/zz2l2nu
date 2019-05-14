# Locates yaml-cpp library
#
# Files of the library are searched for using locations provided in
# environmental variable YAMLCPP_ROOT_DIR and CMake variable yamlcpp_ROOT_DIR,
# as well as standard locations. The target that defines the library is
# exported.

include(CMakeFindDependencyMacro)
find_dependency(Boost)

find_path(yamlcpp_INCLUDE_DIR yaml-cpp/yaml.h
          HINTS $ENV{YAMLCPP_ROOT_DIR}/include ${yamlcpp_ROOT_DIR}/include)

find_library(yamlcpp_LIBRARY NAMES yaml-cpp
             HINTS $ENV{YAMLCPP_ROOT_DIR}/lib ${yamlcpp_ROOT_DIR}/lib)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(yamlcpp
  REQUIRED_VARS yamlcpp_LIBRARY yamlcpp_INCLUDE_DIR)

mark_as_advanced(yamlcpp_FOUND yamlcpp_LIBRARY yamlcpp_INCLUDE_DIR)

if(yamlcpp_FOUND)
  set(yamlcpp_INCLUDE_DIRS ${yamlcpp_INCLUDE_DIR})
  set(yamlcpp_LIBRARIES ${yamlcpp_LIBRARY})
endif()

if(yamlcpp_FOUND AND NOT TARGET yamlcpp::yamlcpp)
  add_library(yamlcpp::yamlcpp SHARED IMPORTED
    ${yamlcpp_LIBRARY}
  )
  set_target_properties(yamlcpp::yamlcpp PROPERTIES
    IMPORTED_LOCATION ${yamlcpp_LIBRARY}
    INTERFACE_INCLUDE_DIRECTORIES ${yamlcpp_INCLUDE_DIR}
    INTERFACE_LINK_LIBRARIES Boost::boost
  )
endif()

