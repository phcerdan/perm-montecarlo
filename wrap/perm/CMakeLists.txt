set(python_module_name _perm)
pybind11_add_module(${python_module_name} MODULE
  perm_common_types_py.cpp
  lattice_lut_py.cpp
  perm_init_py.cpp
  energy_functions_py.cpp
  lattice_boundary_py.cpp
  perm_py.cpp
  )
target_link_libraries(${python_module_name} PRIVATE PERMCore)

install(TARGETS ${python_module_name}
  DESTINATION perm
  COMPONENT runtime
  )

# Copy at build time
add_custom_command(
  TARGET ${python_module_name}
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
  $<TARGET_FILE:${python_module_name}>
  $<TARGET_FILE_NAME:${python_module_name}>
  )

# build tree
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/plot
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
# install tree
install(FILES __init__.py
  DESTINATION perm
  COMPONENT runtime
  )
install(DIRECTORY plot
  DESTINATION perm
  COMPONENT runtime
  )

if(PERM_BUILD_TESTING)
  set(test_folder "${CMAKE_CURRENT_SOURCE_DIR}/test_python")
  set(python_tests_
    test_perm.py
  )

  # test files should start with "test_"
  # unittest functions (in .py) should start with "test_" for discover to work
  foreach(python_test ${python_tests_})
      add_test(NAME python||${python_module_name}||${python_test}
        COMMAND ${PYTHON_EXECUTABLE}
          -m unittest discover
          -s ${test_folder}
          -p ${python_test}
      )
  endforeach()
endif()
