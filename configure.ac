#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
AC_INIT([ransfoil], [2.2], [sdhzhs@gmail.com])
AM_INIT_AUTOMAKE
LT_PREREQ([2.4])
AC_CANONICAL_HOST
AS_CASE([$host],[i?86-pc-mingw32], [LTIOP=win32-dll],[LTIOP=""])
AS_CASE([$host],[x86_64-pc-mingw32], [LTIOP=win64-dll],[LTIOP=""])
LT_INIT([disable-static ${LTIOP}])
AM_CONDITIONAL([BUILDDLL], [test "x${LTIOP}" = xwin32-dll -o "x${LTIOP}" = xwin64-dll])
AC_CONFIG_SRCDIR([src/console.f03])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([config.h])
AC_PROG_FC
AC_PROG_INSTALL
AC_ARG_ENABLE([static],[],[staticlib=${enableval}], [staticlib=no])
AM_CONDITIONAL([BUILDSLL], [test "x${staticlib}" = xyes])
AC_ARG_ENABLE([hyprelib],[AS_HELP_STRING([--enable-hyprelib],[enable hypre library])],[hyprelib=${enableval}], [hyprelib=no])
AM_CONDITIONAL([LINKHYPRE], [test "x${hyprelib}" = xyes])
AS_IF([test "x${hyprelib}" = xyes],[AC_DEFINE([HYPRE],[1],[Define to 1 if you have the HYPRE library])])
AC_ARG_ENABLE([hyprecudalib],[AS_HELP_STRING([--enable-hyprecudalib],[enable hypre library with GPU CUDA support])],[hyprecudalib=${enableval}], [hyprecudalib=no])
AM_CONDITIONAL([LINKHYPRECUDA], [test "x${hyprecudalib}" = xyes])
AS_IF([test "x${hyprecudalib}" = xyes],[AC_DEFINE([HYPRE_GPU],[1],[Define to 1 if you have the HYPRE library with GPU support])])
AC_ARG_WITH([hyprelib-dir],[AS_HELP_STRING([--with-hyprelib-dir=DIR],[DIR is the top directory containing the hypre library if --enable-hyprelib is set])],\
[hypredir=${withval}], [])
AC_SUBST([hypredir])
AC_ARG_ENABLE([openmp],[AS_HELP_STRING([--enable-openmp],[enable openmp parallelization])],[openmp=${enableval}], [openmp=no])
AM_CONDITIONAL([LINKOMP], [test "x${openmp}" = xyes])
AC_ARG_ENABLE([debug], [AS_HELP_STRING([--enable-debug],[enable debug option])], [debug=${enableval}], [debug=no])
AM_CONDITIONAL([DEBUG], [test "x${debug}" = xyes])
AC_ARG_ENABLE([profile], [AS_HELP_STRING([--enable-profile],[enable profile option])], [profile=${enableval}], [profile=no])
AM_CONDITIONAL([PROF], [test "x${profile}" = xyes])
AC_CONFIG_FILES([Makefile lib/Makefile src/Makefile cases/Makefile doc/Makefile scripts/Makefile scripts/subsonic/Makefile scripts/supersonic/Makefile \
scripts/transonic/Makefile scripts/whitcomb1/Makefile])
AC_OUTPUT
