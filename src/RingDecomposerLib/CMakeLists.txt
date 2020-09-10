add_library(
  RingDecomposerLib
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
  RDLbitset.c
)
set_property(TARGET RingDecomposerLib PROPERTY C_STANDARD 90)
set_target_properties(RingDecomposerLib PROPERTIES POSITION_INDEPENDENT_CODE ON)
if(NOT WIN32)
  target_link_libraries(RingDecomposerLib m)
endif()
install(
  TARGETS RingDecomposerLib
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
)

install(FILES RingDecomposerLib.h DESTINATION include)
