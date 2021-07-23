add_executable(Test Test.c TestDemo.c TestValidate.c)
set_property(TARGET Test PROPERTY C_STANDARD 90)
target_link_libraries(Test RingDecomposerLib)
