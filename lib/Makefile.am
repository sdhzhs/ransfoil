AUTOMAKE_OPTIONS = foreign
FCFLAGS = -Wall -ffree-line-length-none -fdefault-real-8 -fdefault-double-8 -J$(top_srcdir)/lib -O3
lib_LTLIBRARIES = libaero2d.la
libaero2d_la_SOURCES = COM.f03 BlockTDMA.f03 Cell2node.f03 Condiff.f03 Connector.f03 Csfit.f03 Defercorrect.f03 Densityinterp.f03 Derivatives.f03 Genmesh.f03 Hypmesh.f03 \
Initial.f03 Interp1d.f03 Jacobi.f03 Libaero2d.F03 Linearsolver.f03 PCEcoe.f03 Postprocess.f03 Readpara.f03 Results.f03 RMS.f03 Turvis.f03 \
Velcorrect.f03 Wallfunc.f03
libaero2d_la_LDFLAGS = -release 2.2
if BUILDDLL
libaero2d_la_LDFLAGS += -no-undefined
endif
if LINKHYPRE
FCFLAGS += -I$(hypredir)/include
libaero2d_la_SOURCES += Hypresolver.f03
libaero2d_la_LDFLAGS += -L$(hypredir)/lib -lHYPRE
endif
if LINKHYPRECUDA
FCFLAGS += -I$(hypredir)/include
libaero2d_la_SOURCES += cudaf.f03 HypreCUDAsolver.f03 
libaero2d_la_LDFLAGS += -L$(hypredir)/lib -lHYPRE
endif
if LINKOMP
FCFLAGS += -fopenmp
libaero2d_la_LDFLAGS += -lgomp
endif
if DEBUG
FCFLAGS += -g -O0 -fcheck=bounds -ffpe-trap=invalid,zero,overflow
endif
if PROF
FCFLAGS += -g -pg
endif
nodist_include_HEADERS = aero2dcom.mod
