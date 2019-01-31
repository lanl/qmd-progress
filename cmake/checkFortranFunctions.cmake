# Check if Fortran functions exist

macro(check_norm2)
  message(STATUS "Looking for Fortran NORM2 function")
  file(WRITE
    ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/testNorm2.f90
    "
      program testNorm2
      real a(3)
      a = 1
      anorm = norm2(a)
      end program testNorm2
    "
    )
    try_compile(HAVE_NORM2
    ${CMAKE_BINARY_DIR}
    ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/testNorm2.f90
    )
    if(HAVE_NORM2)
      message(STATUS "Looking for Fortran NORM2 - found")
      add_definitions(-DNORM2)
    else()
      message(STATUS "Looking for Fortran ${FUNCTION} - not found")
    endif()
endmacro()
