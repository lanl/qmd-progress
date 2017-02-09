#!/bin/bash

RUN="../../build/mdresponse"  #MDResponse program 

time $RUN input.in | tee out

echo -e "\nEnd of run"
