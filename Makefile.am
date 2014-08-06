
ACLOCAL_AMFLAGS = -I m4
SUBDIRS = src
EXTRA_DIST = Rsrc

dist-hook:
	rm -rf `find $(distdir)/Rsrc -name .svn`
	rm -rf `find $(distdir)/Rsrc -name ._*`
	rm -rf `find $(distdir)/Rsrc -name *.o`
	rm -rf `find $(distdir)/Rsrc -name *.so`
	rm -rf `find $(distdir)/m4 -name ._*`

