find_path(
  MelaCore_INCLUDE_DIR Mela.h
  HINTS "$ENV{MELA_ROOT_DIR}/interface"
)
set(MelaCore_LIBRARY_DIR_HINT "$ENV{MELA_ROOT_DIR}/data/$ENV{SCRAM_ARCH}")
find_library(
  MelaCore_JHUGenMELAMELA_LIBRARY JHUGenMELAMELA
  HINTS "${MelaCore_LIBRARY_DIR_HINT}"
)
find_library(
  MelaCore_collier_LIBRARY collier
  HINTS "${MelaCore_LIBRARY_DIR_HINT}"
  NO_DEFAULT_PATH  # Ignore the library in the LCG environment
)
find_library(
  MelaCore_mcfm_LIBRARY NAMES mcfm_707 mcfm
  HINTS "${MelaCore_LIBRARY_DIR_HINT}"
)

find_path(
  MelaAnalytics_CandidateLOCaster_INCLUDE_DIR MELACandidateRecaster.h
  HINTS "$ENV{MELA_ANALYTICS_ROOT_DIR}/CandidateLOCaster/interface"
)
find_library(
  MelaAnalytics_CandidateLOCaster_LIBRARY MelaAnalyticsCandidateLOCaster
  HINTS "$ENV{MELA_ANALYTICS_ROOT_DIR}/CandidateLOCaster/lib"
)

find_path(
  MelaAnalytics_EventContainer_INCLUDE_DIR MELAEvent.h
  HINTS "$ENV{MELA_ANALYTICS_ROOT_DIR}/EventContainer/interface"
)
find_library(
  MelaAnalytics_EventContainer_LIBRARY MelaAnalyticsEventContainer
  HINTS "$ENV{MELA_ANALYTICS_ROOT_DIR}/EventContainer/lib"
)

find_path(
  MelaAnalytics_GenericMEComputer_INCLUDE_DIR GMECHelperFunctions.h
  HINTS "$ENV{MELA_ANALYTICS_ROOT_DIR}/GenericMEComputer/interface"
)
find_library(
  MelaAnalytics_GenericMEComputer_LIBRARY MelaAnalyticsGenericMEComputer
  HINTS "$ENV{MELA_ANALYTICS_ROOT_DIR}/GenericMEComputer/lib"
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Mela
  REQUIRED_VARS
    MelaCore_INCLUDE_DIR MelaCore_JHUGenMELAMELA_LIBRARY
    MelaCore_collier_LIBRARY MelaCore_mcfm_LIBRARY
    MelaAnalytics_CandidateLOCaster_INCLUDE_DIR
    MelaAnalytics_CandidateLOCaster_LIBRARY
    MelaAnalytics_EventContainer_INCLUDE_DIR
    MelaAnalytics_EventContainer_LIBRARY
    MelaAnalytics_GenericMEComputer_INCLUDE_DIR
    MelaAnalytics_GenericMEComputer_LIBRARY
)

mark_as_advanced(
    Mela_FOUND MelaCore_INCLUDE_DIR MelaCore_JHUGenMELAMELA_LIBRARY
    MelaCore_collier_LIBRARY MelaCore_mcfm_LIBRARY
    MelaAnalytics_CandidateLOCaster_INCLUDE_DIR
    MelaAnalytics_CandidateLOCaster_LIBRARY
    MelaAnalytics_EventContainer_INCLUDE_DIR
    MelaAnalytics_EventContainer_LIBRARY
    MelaAnalytics_GenericMEComputer_INCLUDE_DIR
    MelaAnalytics_GenericMEComputer_LIBRARY
)

if(Mela_FOUND AND NOT TARGET Mela)
  add_library(Mela::MelaCore::Mela SHARED IMPORTED)
  set_target_properties(Mela::MelaCore::Mela PROPERTIES
    IMPORTED_LOCATION "${MelaCore_JHUGenMELAMELA_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${MelaCore_INCLUDE_DIR}"
  )

  add_library(Mela::MelaCore::collier SHARED IMPORTED)
  set_target_properties(Mela::MelaCore::collier PROPERTIES
    IMPORTED_LOCATION "${MelaCore_collier_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${MelaCore_INCLUDE_DIR}"
  )

  add_library(Mela::MelaCore::mcfm SHARED IMPORTED)
  set_target_properties(Mela::MelaCore::mcfm PROPERTIES
    IMPORTED_LOCATION "${MelaCore_mcfm_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${MelaCore_INCLUDE_DIR}"
  )

  add_library(Mela::MelaAnalytics::CandidateLOCaster SHARED IMPORTED)
  set_target_properties(Mela::MelaAnalytics::CandidateLOCaster PROPERTIES
    IMPORTED_LOCATION "${MelaAnalytics_CandidateLOCaster_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${MelaAnalytics_CandidateLOCaster_INCLUDE_DIR}"
  )

  add_library(Mela::MelaAnalytics::EventContainer SHARED IMPORTED)
  set_target_properties(Mela::MelaAnalytics::EventContainer PROPERTIES
    IMPORTED_LOCATION "${MelaAnalytics_EventContainer_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${MelaAnalytics_EventContainer_INCLUDE_DIR}"
  )

  add_library(Mela::MelaAnalytics::GenericMEComputer SHARED IMPORTED)
  set_target_properties(Mela::MelaAnalytics::GenericMEComputer PROPERTIES
    IMPORTED_LOCATION "${MelaAnalytics_GenericMEComputer_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${MelaAnalytics_GenericMEComputer_INCLUDE_DIR}"
  )

  add_library(Mela INTERFACE)
  target_link_libraries(Mela
    INTERFACE
      Mela::MelaCore::Mela Mela::MelaCore::collier Mela::MelaCore::mcfm
      Mela::MelaAnalytics::CandidateLOCaster Mela::MelaAnalytics::EventContainer
      Mela::MelaAnalytics::GenericMEComputer
  )
endif()
