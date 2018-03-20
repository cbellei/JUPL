#!/bin/bash

# gfortran -w constants.f90 multION.f90 sources.f90 ibc.f90 -o mION

if ls *.csv 1> /dev/null 2>&1; then
   echo 	
   echo "simulation output files found - do you want to delete them? (y/n)"
   read del
   echo 
   if [ $del == "y" ]; then
      echo "***** deleting files *****"
      rm -f *.csv
   else
      echo "***** aborting simulation *****"
      echo
      exit
   fi
else
   rm -f *.csv
fi

echo
echo "***** starting simulation *****"
echo

#gfortran -fcheck=all -g -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow -Wall -Wtabs constants.f90 multION.f90 io.f90 utils.f90 ibc.f90 -o multION
#gfortran -Wall -Wtabs -Wextra -O3 constants.f90 multION.f90 io.f90 utils.f90 ibc.f90 -o multION


gfortran -fdefault-real-8 -O3 constants.f90 multION.f90 io.f90 utils.f90 ibc.f90 -o multION -llapack -lblas
# gfortran -Wall -Wtabs -Wextra -O3 constants.f90 multION.f90 io.f90 utils.f90 ibc.f90 -o multION -llapack -lblas

if [ $? -ne 0 ]; then
  echo "Errors compiling code"
  exit
fi
./multION
