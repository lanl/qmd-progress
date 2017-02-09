#!/bin/bash

RUN="../../build/mdresponse"  #MDResponse program 

time $RUN input.in > out

echo -e "\nEnd of run"
