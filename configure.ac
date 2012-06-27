#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT([CosmoLib], [1.0], [margie.petkova@gmail.com])
AM_INIT_AUTOMAKE([1.10 no-define foreign])
AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([pow sqrt])

# Get the NR include path
AC_ARG_WITH(nr,AS_HELP_STRING([--with-nr]
        ,[location of nr])
        ,[NR_DIR="$withval"]
        ,[NR_DIR=""])
AC_SUBST(NR_DIR)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
