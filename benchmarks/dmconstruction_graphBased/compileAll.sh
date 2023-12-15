cd /home/cnegre/bml/
./example_build.sh
cd build ; make -j ; make install 
cd /home/cnegre/qmd_progress/
./example_build_correct.sh
cd build ; make -j ; make install 
cd /home/cnegre/qmd-progress/benchmarks/dmconstruction_graphBased/
rm -rf build ; ./build_cmake.sh 
