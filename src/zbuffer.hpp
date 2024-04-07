#ifndef _tgl_zbuffer_h_
#define _tgl_zbuffer_h_

/*
 * Z buffer
 */

#include "zfeatures.h" 

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
  



#define ZB_Z_BITS 16

#define ZB_POINT_Z_FRAC_BITS 14

#define ZB_POINT_S_MIN ( (1<<13) )
#define ZB_POINT_S_MAX ( (1<<22)-(1<<13) )
#define ZB_POINT_T_MIN ( (1<<21) )
#define ZB_POINT_T_MAX ( (1<<30)-(1<<21) )

#define ZB_POINT_RED_MIN ( (1<<10) )
#define ZB_POINT_RED_MAX ( (1<<16)-(1<<10) )
#define ZB_POINT_GREEN_MIN ( (1<<9) )
#define ZB_POINT_GREEN_MAX ( (1<<16)-(1<<9) )
#define ZB_POINT_BLUE_MIN ( (1<<10) )
#define ZB_POINT_BLUE_MAX ( (1<<16)-(1<<10) )
 
/* 16 bit mode */
#define RGB_TO_PIXEL(r,g,b) \
  (((r) & 0xF800) | (((g) >> 5) & 0x07E0) | ((b) >> 11))
typedef uint16_t PIXEL;
#define PSZB 2 
#define PSZSH 4 


/* 32 bpp copy */

#define RGB16_TO_RGB32(p0, p1, v)                     \
  {                                                   \
    uint32_t g, b, gb;                            \
    g = (v & 0x07E007E0) << 5;                        \
    b = (v & 0x001F001F) << 3;                        \
    gb = g | b;                                       \
    p0 = (gb & 0x0000FFFF) | ((v & 0x0000F800) << 8); \
    p1 = (gb >> 16) | ((v & 0xF8000000) >> 8);        \
  }
 

/* 24 bit packed pixel handling */

/* order: RGBR GBRG BRGB */

/* XXX: packed pixel 24 bit support not tested */
/* XXX: big endian case not optimised */

#if BYTE_ORDER == BIG_ENDIAN

#define RGB16_TO_RGB24(p0, p1, p2, v1, v2)                         \
  {                                                                \
    uint32_t r1, g1, b1, gb1, g2, b2, gb2;                     \
    v1 = (v1 << 16) | (v1 >> 16);                                  \
    v2 = (v2 << 16) | (v2 >> 16);                                  \
    r1 = (v1 & 0xF800F800);                                        \
    g1 = (v1 & 0x07E007E0) << 5;                                   \
    b1 = (v1 & 0x001F001F) << 3;                                   \
    gb1 = g1 | b1;                                                 \
    p0 = ((gb1 & 0x0000FFFF) << 8) | (r1 << 16) | (r1 >> 24);      \
    g2 = (v2 & 0x07E007E0) << 5;                                   \
    b2 = (v2 & 0x001F001F) << 3;                                   \
    gb2 = g2 | b2;                                                 \
    p1 = (gb1 & 0xFFFF0000) | (v2 & 0xF800) | ((gb2 >> 8) & 0xff); \
    p2 = (gb2 << 24) | ((v2 & 0xF8000000) >> 8) | (gb2 >> 16);     \
  }

#else

#define RGB16_TO_RGB24(p0, p1, p2, v1, v2)                         \
  {                                                                \
    uint32_t r1, g1, b1, gb1, g2, b2, gb2;                     \
    r1 = (v1 & 0xF800F800);                                        \
    g1 = (v1 & 0x07E007E0) << 5;                                   \
    b1 = (v1 & 0x001F001F) << 3;                                   \
    gb1 = g1 | b1;                                                 \
    p0 = ((gb1 & 0x0000FFFF) << 8) | (r1 << 16) | (r1 >> 24);      \
    g2 = (v2 & 0x07E007E0) << 5;                                   \
    b2 = (v2 & 0x001F001F) << 3;                                   \
    gb2 = g2 | b2;                                                 \
    p1 = (gb1 & 0xFFFF0000) | (v2 & 0xF800) | ((gb2 >> 8) & 0xff); \
    p2 = (gb2 << 24) | ((v2 & 0xF8000000) >> 8) | (gb2 >> 16);     \
  }

#endif


typedef struct {
  int x,y,z;     /* integer coordinates in the zbuffer */
  int s,t;       /* coordinates for the mapping */
  int r,g,b;     /* color indexes */
  
  float sz,tz;   /* temporary coordinates for mapping */
} ZBufferPoint;

 
struct  ZBuffer {
    int xsize,ysize;
    /** PRECISAMOS SE LIVRAR DESSA VARIAVEL CONFUSA*/
    int linesize; /* line size, in bytes */
    int mode;
    
    uint16_t *zbuf;
    PIXEL *pbuf;
    int frame_buffer_allocated;
    
    int nb_colors;
    uint8_t *dctable;
    int *ctable;
    PIXEL *current_texture;

    static inline ZBuffer* open(int xsize, int ysize, int nb_colors,
                 uint8_t* color_indexes, int* color_table,
                 void* frame_buffer) {
      ZBuffer* zb;
      int size;

      zb = (ZBuffer*)gl_malloc(sizeof(ZBuffer));
      if (zb == NULL) return NULL;

      zb->xsize = xsize;
      zb->ysize = ysize;
      zb->linesize = (xsize * PSZB + 3) & ~3;

      size = zb->xsize * zb->ysize * sizeof(uint16_t);

      zb->zbuf = (uint16_t*)gl_malloc(size);
      if (zb->zbuf == NULL)  { 
        gl_free(zb);
        return NULL;
      }
      if (frame_buffer == NULL) {
        zb->pbuf = (PIXEL*)gl_malloc(zb->ysize * zb->linesize);
        if (zb->pbuf == NULL) {
          gl_free(zb->zbuf);
          if (zb->zbuf == NULL)  { 
            gl_free(zb);
            return NULL;
          }
        }
        zb->frame_buffer_allocated = 1;
      } else {
        zb->frame_buffer_allocated = 0;
        zb->pbuf = (PIXEL*)frame_buffer;
      }
    
      zb->current_texture = NULL;

      return zb; 
    }

  static inline void close(ZBuffer* zb) {
    if (zb->frame_buffer_allocated) gl_free(zb->pbuf);

    gl_free(zb->zbuf);
    gl_free(zb);
  }

  inline void clear(bool clear_z, uint16_t z, bool clear_color, uint8_t r, uint8_t g, uint8_t b) {
    PIXEL color;
    uint16_t x,y; 
    PIXEL* pp;
    color = RGB_TO_PIXEL(r, g, b); 
    if (clear_z) {
      memset(zbuf, z, xsize * ysize * sizeof(uint16_t));
    }
    if (clear_color) { 
      for (y = 0; y < ysize; y++) {
        pp=pbuf+ (y*xsize);
        for (x = 0; x < xsize; x++) {
          pp[x]=color;
        }
      }
    }
  }
  
  inline void resize(void* frame_buffer, int xsize, int ysize) {
    int size;

    /* xsize must be a multiple of 4 */
    xsize = xsize & ~3;

    xsize = xsize;
    ysize = ysize;
    linesize = (xsize * PSZB + 3) & ~3;

    size = xsize * ysize * sizeof(uint16_t);

    gl_free(zbuf);
    zbuf = (uint16_t*)gl_malloc(size);

    if (frame_buffer_allocated) gl_free(pbuf);

    if (frame_buffer == NULL) {
      pbuf = (PIXEL*)gl_malloc(ysize * linesize);
      frame_buffer_allocated = 1;
    } else {
      pbuf = (PIXEL*)frame_buffer;
      frame_buffer_allocated = 0;
    }
  }
  inline void copyFrameBuffer(void* buf, int linesize) {
    uint16_t* q;
    uint32_t *p, *p1, v, w0, w1;
    int y, n;

    q = pbuf;
    p1 = (uint32_t*)buf;

    for (y = 0; y < ysize; y++) {
      p = p1;
      n = xsize >> 2;
      do {
        v = *(uint32_t*)q;
        #if BYTE_ORDER == BIG_ENDIAN
              RGB16_TO_RGB32(w1, w0, v);
        #else
              RGB16_TO_RGB32(w0, w1, v);
        #endif
        p[0] = w0;
        p[1] = w1;

        v = *(uint32_t*)(q + 2);
        #if BYTE_ORDER == BIG_ENDIAN
              RGB16_TO_RGB32(w1, w0, v);
        #else
              RGB16_TO_RGB32(w0, w1, v);
        #endif
        p[2] = w0;
        p[3] = w1;

        q += 4;
        p += 4;
      } while (--n > 0);

      p1 = (uint32_t*)((char*)p1 + linesize);
    }
  }
};
  
 
/* zline.c */

void ZB_plot(ZBuffer *zb,ZBufferPoint *p);
void ZB_line(ZBuffer *zb,ZBufferPoint *p1,ZBufferPoint *p2);
void ZB_line_z(ZBuffer * zb, ZBufferPoint * p1, ZBufferPoint * p2);

/* ztriangle.c */

void ZB_setTexture(ZBuffer *zb, PIXEL *texture);

void ZB_fillTriangleFlat(ZBuffer *zb,
		 ZBufferPoint *p1,ZBufferPoint *p2,ZBufferPoint *p3);

void ZB_fillTriangleSmooth(ZBuffer *zb,
		   ZBufferPoint *p1,ZBufferPoint *p2,ZBufferPoint *p3);

void ZB_fillTriangleMapping(ZBuffer *zb,
		    ZBufferPoint *p1,ZBufferPoint *p2,ZBufferPoint *p3);

void ZB_fillTriangleMappingPerspective(ZBuffer *zb,
                    ZBufferPoint *p0,ZBufferPoint *p1,ZBufferPoint *p2);


typedef void (*ZB_fillTriangleFunc)(ZBuffer  *,
	    ZBufferPoint *,ZBufferPoint *,ZBufferPoint *);

 
#endif /* _tgl_zbuffer_h_ */
