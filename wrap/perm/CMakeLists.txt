set(python_module_name _perm)
pybind11_add_module(${python_module_name} MODULE
  perm_init_py.cpp
  perm_py.cpp
  perm_common_types_py.cpp
  )
target_link_libraries(${python_module_name} PRIVATE PERMCore)

install(TARGETS ${python_module_name} DESTINATION perm)

# build tree
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
# install tree
install(FILES __init__.py DESTINATION perm)

# # PythonExtensions is provided by scikit-build: pip install scikit-build
# find_package(PythonExtensions REQUIRED MODULE)
# if(DEFINED PythonExtensions_DIR OR PythonExtensions_DIR)
#   message(STATUS "scikit-build cmake package found in ${PythonExtensions_DIR}")
#   python_extension_module(${python_module_name} MODULE)
# else()
#   message(AUTHOR_WARNING "'scikit-build' not found, python wrappings won't be packaged")
#   message(AUTHOR_WARNING "CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")
# endif()

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
