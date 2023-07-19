#!/bin/bash
export PROJECT_NAME=${CMAKE_BUILD_TYPE:="PROGRESS"}
export PROJECT_VERSION=${PROJECT_VERSION:="master"}
export PROJECT_DESCRIPTION=${PROJECT_DESCRIPTION:="This library is focused on the development of general solvers that are commonly used in quantum chemistry packages."}
sed -i -e 's/@PROJECT_NAME@/'$PROJECT_NAME'/g' Doxyfile.in
sed -i -e 's/@PROJECT_VERSION@/'$PROJECT_VERSION'/g' Doxyfile.in
sed -i -e 's/@PROJECT_DESCRIPTION@/'$PROJECT_DESCRIPTION'/g' Doxyfile.in

mkdir source/_static; doxygen Doxyfile.in; cp -r ./bin/html ./source/_static/doxy 

                                                                                                                                                                                              
                                    
