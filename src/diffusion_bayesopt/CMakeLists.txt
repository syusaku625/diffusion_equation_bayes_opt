set(main_target "Main")

set(compile_options_tmp ${base_compile_options})
string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

add_executable(${main_target} main.cpp)
target_sources(${main_target}
  PRIVATE
    main.cpp
)
target_compile_options(${main_target} PUBLIC ${compile_options_tmp})

target_include_directories(${main_target}
  PUBLIC ${PROJECT_SOURCE_DIR}/include/bayesopt_include
  PUBLIC ${PROJECT_SOURCE_DIR}/include/diffusion_include
)

target_link_libraries(${main_target} PRIVATE
  ${bayesopt_LIBRARY}
  ${Nlopt_LIBRARY}
  ${TP_LIBRARY}
  ${HDF5_LIBRARIES}
  Base
  iomp5
)

install(TARGETS ${main_target} DESTINATION ${CMAKE_INSTALL_PREFIX})
