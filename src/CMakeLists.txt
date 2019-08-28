set(PERMCore_sources
  perm.cpp
  perm_common_types.cpp
  )

add_library(PERMCore ${PERMCore_sources} )
target_include_directories(PERMCore PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  )
set_property(TARGET PERMCore PROPERTY POSITION_INDEPENDENT_CODE ON)
# install targets
install(TARGETS PERMCore
        EXPORT PERMTargets
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )

install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/perm
  PATTERN "*.txt" EXCLUDE)

# export to the build tree
export( TARGETS PERMCore
  NAMESPACE PERM::
  APPEND FILE ${perm_export_file})
