#! /bin/sh
gfortran -Wall -O3 libcallexpf.f03 -I../lib -L../lib/.libs -laero2d -o caller.exe