#!/bin/bash

TOP_DIR="${PWD}"
BUILD_DIR="${BUILD_DIR:=${TOP_DIR}/build}"
INSTALL_DIR="${INSTALL_DIR:=${TOP_DIR}/install}"
LOG_FILE="${TOP_DIR}/build.log"
VERBOSE_MAKEFILE="${VERBOSE_MAKEFILE:=no}"

help() {
    cat <<EOF
Usage:

This script can be used to build and test the progress library.  The script has to
be given a command. Known commands are:

create      - Create the build and install directories ('build' and 'install')
configure   - Configure the build system
compile     - Compile the sources
install     - Install the compiled sources
testing     - Run the test suite
docs        - Generate the API documentation
indent      - Indent the sources
dist        - Generate a tar file (this only works with git)

The following environment variables can be set to influence the configuration
step and the build:

EOF
    set_defaults
    echo "CMAKE_BUILD_TYPE   {Release,Debug}          (default is ${CMAKE_BUILD_TYPE})"
    echo "CMAKE_PREFIX_PATH  Path to include          (default is '${CMAKE_PREFIX_PATH}')"
    echo "CC                 Path to C compiler       (default is ${CC})"
    echo "CXX                Path to C++ compiler     (default is ${CXX})"
    echo "FC                 Path to Fortran compiler (default is ${FC})"
    echo "BML_OPENMP         {yes,no}                 (default is ${BML_OPENMP})"
    echo "PROGRESS_OPENMP    {yes,no}                 (default is ${PROGRESS_OPENMP})"
    echo "PROGRESS_MPI       {yes,no}                 (default is ${PROGRESS_MPI})"
    echo "PROGRESS_TESTING   {yes,no}                 (default is ${PROGRESS_TESTING})"
    echo "PROGRESS_EXAMPLES  {yes,no}                 (default is ${PROGRESS_EXAMPLES})"
    echo "PROGRESS_GRAPHLIB  {yes,no}                 (default is ${PROGRESS_GRAPHLIB})"
    echo "BUILD_DIR          Path to build dir        (default is ${BUILD_DIR})"
    echo "INSTALL_DIR        Path to install dir      (default is ${INSTALL_DIR})"
    echo "EXTRA_FCFLAGS      Extra fortran flags      (default is '${EXTRA_FCFLAGS}')"
    echo "EXTRA_LINK_FLAGS   Any extra link flag      (default is '${EXTRA_LINK_FLAGS}')"
    echo "SANITY_CHECK       Add sanity checks        (default is ${SANITY_CHECK})"
}

set_defaults() {
    : ${CMAKE_BUILD_TYPE:=Release}
    : ${CMAKE_PREFIX_PATH:=""}
    : "${CC:=gcc}"
    : "${CXX:=g++}"
    : "${FC:=gfortran}"
    : ${BML_OPENMP:=yes}
    : ${PROGRESS_OPENMP:=yes}
    : ${PROGRESS_MPI:=no}
    : ${PROGRESS_TESTING:=no}
    : ${PROGRESS_EXAMPLES:=no}
    : ${PROGRESS_GRAPHLIB:=no}
    : "${EXTRA_FCFLAGS:=}"
    : ${EXTRA_LINK_FLAGS:=""}
    : ${SANITY_CHECK:=no}
}

die() {
    echo "fatal error"
    if [[ -f "${BUILD_DIR}/CMakeFiles/CMakeOutput.log" ]]; then
        echo "appending CMake output"
        echo "*********** CMake Output ***********" >> ${LOG_FILE}
        cat "${BUILD_DIR}/CMakeFiles/CMakeOutput.log" >> ${LOG_FILE}
    fi
    if [[ -f "${BUILD_DIR}/CMakeFiles/CMakeError.log" ]]; then
        echo "appending CMake error"
        echo "*********** CMake Error ***********" >> ${LOG_FILE}
        cat "${BUILD_DIR}/CMakeFiles/CMakeError.log" >> ${LOG_FILE}
    fi
    echo "the output from this build was written to ${LOG_FILE}"
    if [[ $# -gt 1 ]]; then
        exit $1
    else
        exit 1
    fi
}

check_pipe_error() {
    if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
        die ${PIPESTATUS[0]}
    fi
}

create() {
    mkdir -v -p "${BUILD_DIR}" || die
    mkdir -v -p "${INSTALL_DIR}" || die
}

configure() {
    set_defaults
    cd "${BUILD_DIR}"
    if [[ -f "${BUILD_DIR}/CMakeCache.txt" ]]; then
        rm -v "${BUILD_DIR}/CMakeCache.txt" || die
    fi
    ${CMAKE:=cmake} .. \
        -DCMAKE_BUILD_TYPE="${CMAKE_BUILD_TYPE}" \
	${CMAKE_PREFIX_PATH:+-DCMAKE_PREFIX_PATH="${CMAKE_PREFIX_PATH}"} \
        -DCMAKE_C_COMPILER="${CC}" \
        -DCMAKE_CXX_COMPILER="${CXX}" \
        -DCMAKE_Fortran_COMPILER="${FC}" \
        ${CMAKE_C_FLAGS:+-DCMAKE_C_FLAGS="${CMAKE_C_FLAGS}"} \
        ${CMAKE_CXX_FLAGS:+-DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS}"} \
        ${CMAKE_Fortran_FLAGS:+-DCMAKE_Fortran_FLAGS="${CMAKE_Fortran_FLAGS}"} \
        -DCMAKE_INSTALL_PREFIX="${INSTALL_DIR}" \
        -DBML_OPENMP="${BML_OPENMP}" \
        -DPROGRESS_OPENMP="${PROGRESS_OPENMP}" \
        -DPROGRESS_MPI="${PROGRESS_MPI}" \
        -DBUILD_SHARED_LIBS="${BUILD_SHARED_LIBS:=no}" \
        -DPROGRESS_TESTING="${PROGRESS_TESTING}" \
        -DPROGRESS_EXAMPLES="${PROGRESS_EXAMPLES}" \
        -DPROGRESS_GRAPHLIB="${PROGRESS_GRAPHLIB}" \
        -DEXTRA_FCFLAGS="${EXTRA_FCFLAGS}" \
        -DEXTRA_LINK_FLAGS="${EXTRA_LINK_FLAGS}" \
        -DCMAKE_VERBOSE_MAKEFILE=${VERBOSE_MAKEFILE} \
        -DSANITY_CHECK=${SANITY_CHECK} \
        | tee -a "${LOG_FILE}"
    check_pipe_error
    cd "${TOP_DIR}"
}

compile() {
    make -C "${BUILD_DIR}" | tee -a "${LOG_FILE}"
    check_pipe_error
}

docs() {
    make -C "${BUILD_DIR}" docs 2>&1 | tee -a "${LOG_FILE}"
    check_pipe_error
    #make -C "${BUILD_DIR}/doc/latex" 2>&1 | tee -a "${LOG_FILE}"
    #check_pipe_error
    #if test -f "${BUILD_DIR}/doc/latex/refman.pdf"; then
    #cp -v "${BUILD_DIR}/doc/latex/refman.pdf" "${TOP_DIR}/bml-manual.pdf"
    #fi
}

install() {
    make -C "${BUILD_DIR}" install 2>&1 | tee -a "${LOG_FILE}"
    check_pipe_error
}

testing() {
    cd "${BUILD_DIR}"
    ctest --output-on-failure 2>&1 | tee -a "${LOG_FILE}"
    check_pipe_error
    cd "${TOP_DIR}"
}

indent() {
    cd "${BUILD_DIR}"
    "${TOP_DIR}/indent.sh" 2>&1 | tee -a "${LOG_FILE}"
    check_pipe_error
    cd "${TOP_DIR}"
    git diff 2>&1 | tee -a "${LOG_FILE}"
    check_pipe_error
    LINES=$(git diff | wc -l)
    if test ${LINES} -gt 0; then
        echo "sources were not formatted correctly"
        die
    fi
}

dist() {
    make -C "${BUILD_DIR}" dist 2>&1 | tee -a "${LOG_FILE}"
    check_pipe_error
}

echo "Writing output to ${LOG_FILE}"

if [[ $# -gt 0 ]]; then
    if [[ "$1" = "-h" || "$1" = "--help" ]]; then
        help
        exit 0
    fi

    case "$1" in
        "create")
            create
            ;;
        "configure")
            create
            configure
            ;;
        "docs")
            create
            configure
            docs
            ;;
        "compile")
            create
            configure
            compile
            ;;
        "install")
            create
            configure
            install
            ;;
        "testing")
            create
            configure
            install
            testing
            ;;
        "indent")
            create
            indent
            ;;
        "dist")
            create
            configure
            dist
            ;;
        *)
            echo "unknown command $1"
            exit 1
            ;;
    esac
else
    echo "missing action"
    help
fi

echo "The output was written to ${LOG_FILE}"
