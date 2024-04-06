include ../config.mk

OBJS= api.o list.o vertex.o init.o matrix.o texture.o \
      misc.o clear.o light.o clip.o select.o get.o error.o \
      zbuffer.o zline.o zdither.o ztriangle.o \
      image_util.o oscontext.o msghandling.o \
      arrays.o specbuf.o memory.o
ifdef TINYGL_USE_GLX
OBJS += glx.o
endif
ifdef TINYGL_USE_NANOX
OBJS += nglx.o
endif

INCLUDES = -I../include
LIB = libTinyGL.a

all: $(LIB)

$(LIB): $(OBJS)
	rm -f $(LIB)
	ar rcs $(LIB) $(OBJS)
	cp $(LIB) ../lib

clean:
	rm -f *~ *.o *.a

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $*.c

clip.o: zgl.hpp zfeatures.h
vertex.o: zgl.hpp zfeatures.h
light.o: zgl.hpp zfeatures.h
matrix.o: zgl.hpp zfeatures.h
list.o: zgl.hpp opinfo.h zfeatures.h
arrays.c: zgl.hpp zfeatures.h
specbuf.o: zgl.hpp zfeatures.h
glx.o: zgl.hpp zfeatures.h
nglx.o: zgl.hpp zfeatures.h
zline.o: zgl.hpp zfeatures.h zline.h

ztriangle.o: ztriangle.c ztriangle.h zgl.hpp zfeatures.h
	$(CC) $(CFLAGS) -Wno-uninitialized $(INCLUDES) -c $*.c
