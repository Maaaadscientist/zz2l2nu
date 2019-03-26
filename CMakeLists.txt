cmake_minimum_required(VERSION 3.11 FATAL_ERROR)
project(hzz2l2nu CXX)

# Require C++17 in all targets
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# External dependencies
find_package(Boost 1.69 REQUIRED COMPONENTS program_options)
find_package(ROOT 6 REQUIRED)


# Convenience library
add_library(hzz2l2nu STATIC
  src/BTagCalibrationStandalone.cc
  src/BTagWeight.cc
  src/ElectronBuilder.cc
  src/EWCorrectionWeight.cc
  src/LooperMain.cc
  src/MuonBuilder.cc
  src/NRBLooper.cc
  src/ObjectSelection.cc
  src/Options.cc
  src/PhotonBuilder.cc
  src/PileUpWeight.cc
  src/RoccoR.cc
  src/SmartSelectionMonitor.cc
  src/SmartSelectionMonitor_hzz.cc
  src/Tables.cc
  src/Trigger.cc
  src/Utils.cc
)
target_include_directories(hzz2l2nu PUBLIC include)
target_link_libraries(hzz2l2nu
  PUBLIC Boost::boost Boost::program_options
  PUBLIC ROOT::Hist ROOT::MathCore ROOT::Physics ROOT::Tree ROOT::TreePlayer
)

add_executable(runHZZanalysis src/runHZZanalysis.cc)
target_link_libraries(runHZZanalysis PRIVATE hzz2l2nu Boost::boost)


# Will install the executable in the source directory
set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}")
install(TARGETS runHZZanalysis DESTINATION bin)

