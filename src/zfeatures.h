#ifndef _tgl_features_h_
#define _tgl_features_h_

#include "msghandling.hpp"

/* It is possible to enable/disable (compile time) features in this
   header file. */

#define TGL_FEATURE_ARRAYS         1
#define TGL_FEATURE_DISPLAYLISTS   1
#define TGL_FEATURE_POLYGON_OFFSET 1

/*
 * Matrix of internal and external pixel formats supported. 'Y' means
 * supported.
 * 
 *           External  8    16    24    32
 * Internal 
 *  16                 Y     Y     Y     Y
 * 
 *
 * 15 bpp does not work yet (although it is easy to add it - ask me if
 * you need it).
 * 
 * Internal pixel format: 16BITS
 * External pixel format: see TGL_FEATURE_xxx_BITS 
 */

/* enable various convertion code from internal pixel format (usually
   16 bits per pixel) to any external format */
#define TGL_FEATURE_16_BITS        1
#define TGL_FEATURE_8_BITS         1
#define TGL_FEATURE_24_BITS        1
#define TGL_FEATURE_32_BITS        1

/*
 * Memory allocator for TinyGL
 */

#include <stdlib.h>
/* modify these functions so that they suit your needs */

static inline void gl_free(void *p) {
    free(p);
}

static inline void *gl_malloc(int size) {
    return malloc(size);
}

static inline void *gl_zalloc(int size) {
    return calloc(1, size);
}

#include <stdarg.h> 

static inline void gl_fatal_error(char *format, ...) {
  va_list ap;

  va_start(ap,format);

  fprintf(stderr,"TinyGL: fatal error: ");
  vfprintf(stderr,format,ap);
  fprintf(stderr,"\n");
  exit(1);

  va_end(ap);
}



#endif /* _tgl_features_h_ */
