# This makefile is for cablat, which must be unzipped in a blatSrc/ build.

include ../inc/common.mk
CXX=g++
%.o: %.cpp
	${CXX} ${COPT} ${CFLAGS} ${HG_DEFS} ${HG_WARN} ${HG_INC} ${XINC} -o $@ -c $<

# better to export this from command line... MACHTYPE also used in ../inc/common.mk
# MACHTYPE = x86_64

L = -lm $(SOCKETLIB)
MYLIBDIR = ../lib/$(MACHTYPE)
MYLIBS =  $(MYLIBDIR)/jkOwnLib.a $(MYLIBDIR)/jkweb.a 

O = cablat.o

PROFILE = #-pg -fno-inline

cablat: $O $(MYLIBS)
#	${CC} ${COPT} ${CFLAGS} -o ${BINDIR}/blat $O $(MYLIBS) $L
	${CXX} ${PROFILE} ${COPT} ${CFLAGS} -o cablat $O $(MYLIBS) $L
	${STRIP} cablat

all:
	cd ../lib && ${MAKE}
	make

tags:
	ctags *.c *.h ../lib/*.c ../inc/*.h

clean::
	rm -f $(O) cablat
