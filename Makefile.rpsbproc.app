#$Id: Makefile.rpsbproc.app Shennan Lu $

APP = rpsbproc
SRC = rpsbproc

CPPFLAGS = -I. \
	$(ORIG_CPPFLAGS)

LIB = blastxml xser xutil xncbi

LIBS = $(DL_LIBS) $(ORIG_LIBS)
