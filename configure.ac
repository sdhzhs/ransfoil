#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
AC_INIT([ransfoil], [2.0], [sdhzhs@gmail.com])
AM_INIT_AUTOMAKE(ransfoil, 2.0, sdhzhs@gmail.com)
LT_PREREQ([2.4])
AC_CANONICAL_HOST
AS_CASE([$host],[i?86-pc-mingw32], [LTIOP=win32-dll],[LTIOP=""])
LT_INIT([disable-static ${LTIOP}])
AM_CONDITIONAL([BUILDDLL], [test "x${LTIOP}" = xwin32-dll])
AC_CONFIG_SRCDIR([console.f03])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([config.h])
AC_PROG_FC
AC_PROG_INSTALL
AC_ARG_ENABLE([hyprelib],[AS_HELP_STRING([--enable-hyprelib],[enable hypre library])],[hyprelib=${enableval}], [hyprelib=no])
AM_CONDITIONAL([LINKHYPRE], [test "x${hyprelib}" = xyes])
AS_IF([test "x${hyprelib}" = xyes],[AC_DEFINE([HYPRE],[1],[Define to 1 if you have the HYPRE library])])
AC_ARG_WITH([hyprelib-dir],[AS_HELP_STRING([--with-hyprelib-dir=DIR],[DIR is the top directory containing the hypre library if --enable-hyprelib is set])],\
[hypredir=${withval}], [])
AC_SUBST([hypredir])
AC_ARG_ENABLE([openmp],[AS_HELP_STRING([--enable-openmp],[enable openmp parallelization])],[openmp=${enableval}], [openmp=no])
AM_CONDITIONAL([LINKOMP], [test "x${openmp}" = xyes])
AC_CONFIG_FILES([Makefile case1/Makefile doc/Makefile lib/Makefile])
AC_OUTPUT