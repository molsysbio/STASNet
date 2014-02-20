
ACLOCAL_AMFLAGS = -I m4
SUBDIRS = src
EXTRA_DIST = r_binding

dist-hook:
	rm -rf `find $(distdir)/r_binding -name .svn`
	rm -rf `find $(distdir)/r_binding -name ._*`
	rm -rf `find $(distdir)/r_binding -name *.o`
	rm -rf `find $(distdir)/r_binding -name *.so`
	rm -rf `find $(distdir)/m4 -name ._*`


