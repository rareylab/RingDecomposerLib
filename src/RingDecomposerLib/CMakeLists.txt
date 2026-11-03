function(build_library LIBRARY_NAME SPECIFIER)
  add_library("${LIBRARY_NAME}" "${SPECIFIER}"
              RDLutility.c
              RDLhandler.c
              RDLrelation.c
              RDLapsp.c
              RingDecomposerLib.c
              RDLgraph.c
              RDLinfo.c
              RDLcycleFams.c
              RDLstack.c
              RDLtesting.c
              RDLtarjan.c
              RDLdimacs.c
              RDLbitset.c)
  set_property(TARGET "${LIBRARY_NAME}" PROPERTY C_STANDARD 90)
  set_target_properties("${LIBRARY_NAME}" PROPERTIES POSITION_INDEPENDENT_CODE True)
  if(NOT WIN32)
    target_link_libraries("${LIBRARY_NAME}" m)
  endif()
  install(TARGETS ${LIBRARY_NAME}
          ARCHIVE DESTINATION lib
          LIBRARY DESTINATION lib
          RUNTIME DESTINATION bin)
endfunction()

build_library("${SHARED_LIBRARY}" "SHARED")
build_library("${STATIC_LIBRARY}" "STATIC")

install(FILES RingDecomposerLib.h
        DESTINATION include)

option(RDL_BUILD_ALLOCATION_TESTS "Build allocation failure injection tests" OFF)
if(RDL_BUILD_ALLOCATION_TESTS)
  get_target_property(RDL_ALLOCATION_SOURCES ${STATIC_LIBRARY} SOURCES)
  add_library(RDLAllocationTestLibrary STATIC ${RDL_ALLOCATION_SOURCES})
  target_compile_definitions(RDLAllocationTestLibrary PRIVATE
    malloc=audit_malloc calloc=audit_calloc realloc=audit_realloc free=audit_free)
  if(MSVC)
    target_compile_options(RDLAllocationTestLibrary PRIVATE
      /FI${CMAKE_CURRENT_SOURCE_DIR}/tests/failure_alloc.h)
  else()
    target_compile_options(RDLAllocationTestLibrary PRIVATE
      -include ${CMAKE_CURRENT_SOURCE_DIR}/tests/failure_alloc.h)
  endif()
  add_executable(RDLAllocationFailures tests/allocation_failures.c)
  target_include_directories(RDLAllocationFailures PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(RDLAllocationFailures RDLAllocationTestLibrary)
  if(NOT WIN32)
    target_link_libraries(RDLAllocationFailures m)
  endif()
  add_test(NAME allocation_failures COMMAND RDLAllocationFailures
    ${TEST_DIR}/molecule_00_0.dimacs)
endif()
