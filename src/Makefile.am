SUBDIRS = levmar-2.5

lib_LTLIBRARIES = libfitmodel.la

libfitmodel_la_LIBADD = levmar-2.5/liblevmar.la
libfitmodel_la_SOURCES = fitmodel_in_CPP.cpp  identifiability.cpp mathtree.cpp model.cpp helper_types.cpp generate_response.cpp fitmodel_in_CPP.hpp identifiability.hpp mathtree.hpp model.hpp helper_types.hpp generate_response.hpp 
libfitmodel_la_CPPFLAGS = $(GINAC_CFLAGS) $(BOOST_CPPFLAGS)

libfitmodel_la_includedir=$(includedir)/fitmodel
libfitmodel_la_include_HEADERS = fitmodel_in_CPP.hpp helper_types.hpp mathtree.hpp generate_response.hpp identifiability.hpp model.hpp rref.hpp

dist-hook:
	rm -rf `find $(distdir)/threadpool -name .svn`
