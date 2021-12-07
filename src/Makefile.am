AUTOMAKE_OPTIONS = foreign
FCFLAGS = -Wall -ffree-line-length-none -O3
if DEBUG
FCFLAGS += -g -O0 -fcheck=bounds -ffpe-trap=invalid,zero,overflow
endif
if PROF
FCFLAGS += -g -pg
endif
bin_PROGRAMS = ransfoil
ransfoil_SOURCES = console.f03
ransfoil_LDADD = $(top_srcdir)/lib/libaero2d.la
EXTRA_PROGRAMS = libcallexp
libcallexp_SOURCES = libcallexpf.f03 libcallexpc.c
libcallexp_LDADD = $(top_srcdir)/lib/libaero2d.la
dist_include_HEADERS = aero2dcom.h
