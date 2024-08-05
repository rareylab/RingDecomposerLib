function(build_minimalexample MINIMALEXAMPLE_EXECUTABLE LIBRARY)
  add_executable("${MINIMALEXAMPLE_EXECUTABLE}" MinimalExample.c)
  set (EXTRA_LIBS "${LIBRARY}")
  target_link_libraries ("${MINIMALEXAMPLE_EXECUTABLE}" ${EXTRA_LIBS})
endfunction()

build_minimalexample("${SHARED_MINIMALEXAMPLE}" "${SHARED_LIBRARY}")
build_minimalexample("${STATIC_MINIMALEXAMPLE}" "${STATIC_LIBRARY}")
