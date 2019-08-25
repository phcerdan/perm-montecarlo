# Minimum requirement is a list of cpp:
# set(MODULE_${MODULE_NAME}_TESTS
#  test_get_vtk_points_from_graph.cpp
#  test_graph_points_locator.cpp
# )
# And the dependencies
#  set(MODULE_${MODULE_NAME}_TEST_DEPENDS
#    ${MODULE_${MODULE_NAME}_LIBRARY}
#    ${MODULE_${MODULE_NAME}_DEPENDS}
#    ${GTEST_LIBRARIES})
macro(perm_add_gtests)
  foreach(test_file ${MODULE_${MODULE_NAME}_TESTS})
    string(REGEX REPLACE "\\.[^.]*$" "" test_name "${test_file}")
    add_executable(${test_name} ${test_file})
    target_link_libraries(${test_name} PRIVATE
      ${MODULE_${MODULE_NAME}_TEST_DEPENDS}
      )
    target_include_directories(${test_name} SYSTEM PRIVATE
      ${MODULE_${MODULE_NAME}_TEST_SYSTEM_INCLUDE_DIRS})
    gtest_discover_tests(
      ${test_name}
      TEST_PREFIX ${MODULE_NAME}||${test_name}||
      PROPERTIES LABELS ${MODULE_NAME}
      )
  endforeach()
endmacro(perm_add_gtests)
