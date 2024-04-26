# FindMETIS.cmake
find_path(METIS_INCLUDE_DIR NAMES metis.h)
find_library(METIS_LIBRARY NAMES metis)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Metis DEFAULT_MSG METIS_LIBRARY METIS_INCLUDE_DIR)

if(METIS_FOUND AND NOT TARGET Metis::Metis)
  add_library(Metis::Metis UNKNOWN IMPORTED)
  set_target_properties(Metis::Metis PROPERTIES
    IMPORTED_LOCATION "${METIS_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${METIS_INCLUDE_DIR}"
  )
endif()