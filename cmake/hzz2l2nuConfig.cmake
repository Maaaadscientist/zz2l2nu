get_filename_component(HZZ2L2NU_CMAKE_DIR
  "${CMAKE_CURRENT_LIST_FILE}" DIRECTORY
)
list(APPEND CMAKE_MODULE_PATH "${HZZ2L2NU_CMAKE_DIR}")  # Needed to find yamlcpp

# Exposed dependencies
include(CMakeFindDependencyMacro)
set(Boost_NO_BOOST_CMAKE ON)
find_dependency(Boost 1.70 COMPONENTS log program_options stacktrace_basic)
find_dependency(ROOT 6)
find_dependency(yamlcpp)

if (NOT TARGET hzz2l2nu::hzz2l2nu)
  # The file with exported targets below is generated by CMake
  include("${HZZ2L2NU_CMAKE_DIR}/hzz2l2nuTargets.cmake")
endif()
