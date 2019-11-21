cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

# By default, set the build type to Release
if(NOT DEFINED CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(default_build_type Release)
  message("Using default build type ${default_build_type}")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
    STRING "Select build type" FORCE
  )
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY
    STRINGS "" Debug Release MinSizeRel RelWithDebInfo
  )
endif()

project(hzz2l2nu CXX)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")

# Require C++17 in all targets
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Boost. The first command is needed as a work-around for an issue in
# interaction between CMake and Boost [1].
# [1] https://gitlab.kitware.com/cmake/cmake/issues/18865#note_511829
set(Boost_NO_BOOST_CMAKE ON)
find_package(
  Boost 1.70 REQUIRED
  COMPONENTS log program_options stacktrace_basic
)

# Other external dependencies
find_package(ROOT 6 REQUIRED)
find_package(yamlcpp REQUIRED)


# Minimalistic library to embed the hash of the current Git commit into the
# binary code as its version. CMake will be rerun whenever .git/index changes.
set_property(DIRECTORY APPEND
  PROPERTY CMAKE_CONFIGURE_DEPENDS
  .git/index
)

find_package(Git QUIET REQUIRED)
execute_process(
  COMMAND "${GIT_EXECUTABLE}" rev-parse HEAD
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
  OUTPUT_VARIABLE GIT_COMMIT
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
configure_file(src/Version.cc.in "${CMAKE_BINARY_DIR}/Version.cc")

add_library(version STATIC
  "${CMAKE_BINARY_DIR}/Version.cc"
)
target_include_directories(version PRIVATE include)


# Jet calibration library
add_library(jerc STATIC
  src/JERC/FactorizedJetCorrector.cc
  src/JERC/JetCorrectionUncertainty.cc
  src/JERC/JetCorrectorParameters.cc
  src/JERC/JetResolution.cc
  src/JERC/JetResolutionObject.cc
  src/JERC/SimpleJetCorrectionUncertainty.cc
  src/JERC/SimpleJetCorrector.cc
)
target_include_directories(jerc PRIVATE src/JERC)
target_link_libraries(jerc
  PRIVATE ROOT::GenVector ROOT::Hist
)

# Library to access b tagging scale factors
add_library(btag STATIC
  src/BTag/BTagCalibrationStandalone.cc
)
target_include_directories(btag PRIVATE src/BTag)
target_link_libraries(btag
  PRIVATE ROOT::Hist
)

# Convenience library
add_library(hzz2l2nu STATIC
  src/AnalysisCommon.cc
  src/BTagger.cc
  src/BTagWeight.cc
  src/CollectionBuilder.cc
  src/Dataset.cc
  src/DileptonTrees.cc
  src/ElectronBuilder.cc
  src/EWCorrectionWeight.cc
  src/FileInPath.cc
  src/GenJetBuilder.cc
  src/GenWeight.cc
  src/GenZZBuilder.cc
  src/InstrMetAnalysis.cc
  src/JetBuilder.cc
  src/KFactorCorrection.cc
  src/LeptonWeight.cc
  src/Logger.cc
  src/MainAnalysis.cc
  src/MelaWeight.cc
  src/MeKinFilter.cc
  src/MetFilters.cc
  src/MuonBuilder.cc
  src/NrbAnalysis.cc
  src/ObjectSelection.cc
  src/Options.cc
  src/PhotonBuilder.cc
  src/PileUpWeight.cc
  src/PtMissBuilder.cc
  src/RoccoR.cc
  src/SmartSelectionMonitor.cc
  src/SmartSelectionMonitor_hzz.cc
  src/Tables.cc
  src/TabulatedRandomGenerator.cc
  src/Trigger.cc
  src/Utils.cc
)
target_include_directories(hzz2l2nu PUBLIC include)
target_link_libraries(hzz2l2nu
  PRIVATE stdc++fs
  PRIVATE jerc btag
  PRIVATE version
  PUBLIC Boost::boost Boost::log Boost::program_options
  PUBLIC Boost::stacktrace_basic -rdynamic  # To preserve human-readable names
  PUBLIC ROOT::Hist ROOT::MathCore ROOT::Physics ROOT::Tree ROOT::TreePlayer
  PUBLIC yamlcpp::yamlcpp
)

add_executable(runHZZanalysis src/runHZZanalysis.cc)
target_link_libraries(runHZZanalysis PRIVATE hzz2l2nu Boost::boost)

add_executable(haddws src/haddws.cc)
target_link_libraries(haddws PRIVATE ROOT::Core ROOT::Hist ROOT::RIO ROOT::Tree)

