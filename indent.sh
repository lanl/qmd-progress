#!/bin/bash
#
# Install emacs: 
#	  sudo apt-get install emacs
#
# Get the indentation emacs script as follows: 
#  	wget https://raw.github.com/Glavin001/atom-beautify/master/src/beautifiers/fortran-beautifier/emacs-fortran-formating-script.lisp
# 
# Change to downcase: 
#   sed 's/-upcase-/-downcase-/g' emacs-fortran-formating-script.lisp > emacs-fortran-formating-script-down.lisp
#
# To use it just call this script from the main qmd-progress folder. 
# 	./tools/indent.sh

PROGRESS_PATH=$(dirname $(readlink -f $0))

for file in $PROGRESS_PATH/src/*.F90
do 
  emacs -batch -l $PROGRESS_PATH/emacs-fortran-formating-script-down.lisp -f f90-batch-indent-region $file
done	


