#
# DLCoal
# Matt Rasmussen
# Copyright 2012
#
# Makefile
#

# install prefix paths
prefix = /usr


# C++ compiler options
CXX = g++

CFLAGS := $(CFLAGS) \
    -Wall -fPIC \
    -Isrc


#=============================================================================
# optional CFLAGS

# profiling
ifdef PROFILE
	CFLAGS := $(CFLAGS) -pg
endif

# debugging
ifdef DEBUG	
	CFLAGS := $(CFLAGS) -g
else
	CFLAGS := $(CFLAGS) -O3
endif


#=============================================================================
# DLCoal files

# package
PKG_VERSION:=$(shell python -c 'import dlcoal; print dlcoal.PROGRAM_VERSION_TEXT')
PKG_NAME=dlcoal
PKG=dist/$(PKG_NAME)-$(PKG_VERSION).tar.gz
PKG_DIR=dist/$(PKG_NAME)-$(PKG_VERSION)

# program files
SCRIPTS =  bin/dlcoal_recon \
           bin/dlcoal_sim \
           bin/mpr \
           bin/view_recon
BINARIES = $(SCRIPTS)

DLCOAL_SRC = \
    src/coal.cpp \
    src/duploss.cpp \
    src/itree.cpp \
    src/spidir/birthdeath.cpp \
    src/spidir/common.cpp \
    src/spidir/logging.cpp \
    src/spidir/phylogeny.cpp \
    src/spidir/Tree.cpp \
    src/spidir/top_prior.cpp \



DLCOAL_OBJS = $(DLCOAL_SRC:.cpp=.o)

LIBS =
# `gsl-config --libs`
#-lgsl -lgslcblas -lm


#=======================
# DLCoal C-library files
LIBDLCOAL = lib/libdlcoal.a
LIBDLCOAL_SHARED_NAME = libdlcoal.so
LIBDLCOAL_SHARED = lib/$(LIBDLCOAL_SHARED_NAME)
LIBDLCOAL_SHARED_INSTALL = $(prefix)/lib/$(LIBDLCOAL_SHARED_NAME)
LIBDLCOAL_OBJS = $(DLCOAL_OBJS)


#=============================================================================
# targets

# default targets
all: $(LIBDLCOAL) $(LIBDLCOAL_SHARED)


#-----------------------------
# DLCOAL C-library
lib: $(LIBDLCOAL) $(LIBDLCOAL_SHARED)

$(LIBDLCOAL): $(LIBDLCOAL_OBJS)
	mkdir -p lib
	$(AR) -r $(LIBDLCOAL) $(LIBDLCOAL_OBJS)

$(LIBDLCOAL_SHARED): $(LIBDLCOAL_OBJS) 
	mkdir -p lib
	$(CXX) -o $(LIBDLCOAL_SHARED) -shared $(LIBDLCOAL_OBJS) $(LIBS)


#-----------------------------
# packaging

pkg:
	python make-pkg.py $(PKG_DIR)

$(PKG):
	python make-pkg.py $(PKG_DIR)

#-----------------------------
# install

install: $(BINARIES) $(LIBDLCOAL_SHARED_INSTALL)
	mkdir -p $(prefix)/bin
	cp $(BINARIES) $(prefix)/bin
	echo $(LIBDLCOAL_SHARED_INSTALL)
	python setup.py install --prefix=$(prefix)

pylib: $(LIBDLCOAL_SHARED_INSTALL)
	python setup.py install --prefix=$(prefix)


$(LIBDLCOAL_SHARED_INSTALL): $(LIBDLCOAL_SHARED)
	mkdir -p $(prefix)/lib
	cp $(LIBDLCOAL_SHARED) $(LIBDLCOAL_SHARED_INSTALL)

#=============================================================================
# basic rules

$(DLCOAL_OBJS): %.o: %.cpp
	$(CXX) -c $(CFLAGS) -o $@ $<


clean:
	rm -f $(DLCOAL_OBJS) $(LIBDLCOAL) $(LIBDLCOAL_SHARED)

clean-obj:
	rm -f $(DLCOAL_OBJS)


#=============================================================================
# dependencies

dep:
	touch Makefile.dep
	makedepend -f Makefile.dep -Y src/*.cpp src/*.h

Makefile.dep:
	touch Makefile.dep
	makedepend -f Makefile.dep -Y src/*.cpp src/*.h

include Makefile.dep
# DO NOT DELETE
