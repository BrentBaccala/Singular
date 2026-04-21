ACLOCAL_AMFLAGS = -I ../../m4

if SI_BUILTIN_SINGRUST
  noinst_LTLIBRARIES=singrust.la
  P_PROCS_CPPFLAGS_COMMON = -DSTATIC_VERSION
  P_PROCS_MODULE_LDFLAGS  = -module
else
  module_LTLIBRARIES=singrust.la
  moduledir = $(libdir)/singular/MOD
  P_PROCS_CPPFLAGS_COMMON = -DDYNAMIC_VERSION
  P_PROCS_MODULE_LDFLAGS = -module -export-dynamic -avoid-version
endif

MYINCLUDES =  -I${top_srcdir} -I${top_builddir} \
-I${top_srcdir}/libpolys -I${top_builddir}/libpolys \
$(FACTORY_INCLUDES) $(RESOURCES_INCLUDES) $(OMALLOC_INCLUDES) \
$(FLINT_CFLAGS) $(NTL_CXXFLAGS) $(NTL_CPPFLAGS) $(GMP_CPPFLAGS) \
$(RUSTGB_CFLAGS) $(CPPFLAGS)

singrust_la_SOURCES  = singrust.cc
singrust_la_CPPFLAGS = ${MYINCLUDES} ${P_PROCS_CPPFLAGS_COMMON}
singrust_la_LDFLAGS  = ${P_PROCS_MODULE_LDFLAGS} ${RUSTGB_LIBS} $(SINGULAR_LDFLAGS)

if SING_WIN
if !SI_BUILTIN_SINGRUST
 singrust_la_LIBADD = ${top_builddir}/Singular/libSingular.la
endif
endif
