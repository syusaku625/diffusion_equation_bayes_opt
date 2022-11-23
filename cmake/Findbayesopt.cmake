# httts://qiita.com/shohirose/items/d9bda00a39a113965c5c

find_path(bayesopt_INCLUDE_DIR nlopt.h
  PATHS
    /usr
    /usr/local
    ${bayesopt_DIR}
    $ENV{bayesopt_DIR}
  PATH_SUFFIXES
    include
)

find_library(bayesopt_LIBRARY 
NAMES bayesopt
PATHS
    /usr
    /usr/local
    ${bayesopt_DIR}
    $ENV{bayesopt_DIR}
  PATH_SUFFIXES
    lib
) 

find_library(Nlopt_LIBRARY 
NAMES nlopt
PATHS
    /usr
    /usr/local
    ${bayesopt_DIR}
    $ENV{bayesopt_DIR}
  PATH_SUFFIXES
    lib
)

mark_as_advanced(
  bayesopt_INCLUDE_DIR
  bayesopt_LIBRARY     # ヘッダーのみのライブラリの場合は不要
  )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(bayesopt
  REQUIRED_VARS
    bayesopt_INCLUDE_DIR
    bayesopt_LIBRARY      # ヘッダーのみのライブラリの場合は不要
  )

if(bayesopt_FOUND AND NOT TARGET bayesopt)
  add_library(bayesopt UNKNOWN IMPORTED)
  set_target_properties(bayesopt PROPERTIES
    IMPORTED_LINK_INTERFACE_LANGUAGES ["C"|"CXX"]  # ヘッダーのみのライブラリの場合は不要
    IMPORTED_LOCATION "${bayesopt_LIBRARY}"       # ヘッダーのみのライブラリの場合は不要
    INTERFACE_INCLUDE_DIRECTORIES "${bayesopt_INCLUDE_DIR}"
    )
endif()