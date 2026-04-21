dnl -*- autoconf -*-
dnl
dnl LB_CHECK_RUSTGB — probe for the rustgb cdylib / header.
dnl
dnl Layout expected under the --with-rustgb prefix:
dnl   ${prefix}/include/rustgb.h
dnl   ${prefix}/target/release/librustgb.so
dnl
dnl This matches the in-source layout of the rustgb Cargo crate,
dnl which is what we pass in during development:
dnl   --with-rustgb=$HOME/rustgb
dnl
dnl On success:
dnl   * AC_DEFINE(HAVE_RUSTGB, 1, ...)
dnl   * AC_SUBST(RUSTGB_CFLAGS, RUSTGB_LDFLAGS, RUSTGB_LIBS)
dnl   * AM_CONDITIONAL(SI_BUILTIN_SINGRUST)
dnl On failure (still with-rustgb != no): configure error.
dnl
AC_DEFUN([LB_CHECK_RUSTGB],
[
  AC_ARG_WITH(rustgb,
    [AS_HELP_STRING([--with-rustgb=yes|no|prefix],
      [Use the rustgb Groebner basis engine. Default is no. When given a prefix, expects ${prefix}/include/rustgb.h and ${prefix}/target/release/librustgb.so.])],
    [], [with_rustgb="no"])

  AS_IF([test "x$with_rustgb" != "xno"],
  [
    AS_IF([test "x$with_rustgb" != "xyes"],
      [RUSTGB_PREFIX="$with_rustgb"
       RUSTGB_CFLAGS="-I${RUSTGB_PREFIX}/include"
       RUSTGB_LDFLAGS="-L${RUSTGB_PREFIX}/target/release"],
      [RUSTGB_CFLAGS=""
       RUSTGB_LDFLAGS=""])

    BACKUP_LIBS=${LIBS}
    BACKUP_LDFLAGS=${LDFLAGS}
    BACKUP_CPPFLAGS=${CPPFLAGS}

    CPPFLAGS="${CPPFLAGS} ${RUSTGB_CFLAGS}"
    LDFLAGS="${LDFLAGS} ${RUSTGB_LDFLAGS}"
    LIBS="-lrustgb ${LIBS}"

    AC_CHECK_HEADER([rustgb.h], [],
      [AC_MSG_ERROR([Cannot find rustgb.h. Pass --with-rustgb=PREFIX.])])

    AC_MSG_CHECKING([for rustgb_version in -lrustgb])
    AC_LINK_IFELSE(
      [AC_LANG_PROGRAM(
        [[#include <rustgb.h>]],
        [[(void) rustgb_version();]])],
      [AC_MSG_RESULT([yes])],
      [AC_MSG_RESULT([no])
       AC_MSG_ERROR([Cannot link against librustgb. Check --with-rustgb=PREFIX and ensure librustgb.so is built (cargo build --release in the rustgb crate).])])

    LIBS=${BACKUP_LIBS}
    LDFLAGS=${BACKUP_LDFLAGS}
    CPPFLAGS=${BACKUP_CPPFLAGS}

    RUSTGB_LIBS="${RUSTGB_LDFLAGS} -lrustgb"
    AC_SUBST(RUSTGB_CFLAGS)
    AC_SUBST(RUSTGB_LDFLAGS)
    AC_SUBST(RUSTGB_LIBS)
    AC_DEFINE(HAVE_RUSTGB, 1, [Define if rustgb is to be used])
  ])

  AM_CONDITIONAL([ENABLE_SINGRUST_MODULE], [test "x$with_rustgb" != "xno"])
])
