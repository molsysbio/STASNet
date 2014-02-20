AC_INIT([Network_Reverse_Engineering C++], [1.0], [nils@sys-bio.net],
             [network_reverse_engineering], [http://www.sys-bio.net/])
AC_PREREQ([2.59])
AM_INIT_AUTOMAKE([1.0])
# AM_EXTRA_RECURSIVE_TARGETS([src src/levmar-2.5/])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_F77
LT_INIT


# Checks for header files.
AC_CHECK_HEADERS([float.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_CHECK_FUNCS([pow sqrt])

# Checks for libraries.

AX_LAPACK([],[ AC_MSG_ERROR([LAPACK is required to build..]) ])
AX_BLAS
AX_BOOST_BASE([1.4],[],[ AC_MSG_ERROR([BOOST is required to build..]) ])

PKG_CHECK_MODULES(GINAC, ginac >= 1.3.7)


# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

AC_CONFIG_FILES([Makefile src/Makefile src/levmar-2.5/Makefile])

AC_OUTPUT


