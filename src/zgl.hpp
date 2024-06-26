#ifndef _tgl_zgl_h_
#define _tgl_zgl_h_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <GL/gl.h>
#include "zbuffer.hpp"
#include "zmath.hpp"
#include "zfeatures.h"

#define DEBUG
/* #define NDEBUG */

enum {

#define ADD_OP(a,b,c) OP_ ## a ,

#include "opinfo.h"

};

/* initially # of allocated GLVertexes (will grow when necessary) */
#define POLYGON_MAX_VERTEX 16

/* Max # of specular light pow buffers */
#define MAX_SPECULAR_BUFFERS 8
/* # of entries in specular buffer */
#define SPECULAR_BUFFER_SIZE 1024
/* specular buffer granularity */
#define SPECULAR_BUFFER_RESOLUTION 1024


#define MAX_MODELVIEW_STACK_DEPTH  32
#define MAX_PROJECTION_STACK_DEPTH 8
#define MAX_TEXTURE_STACK_DEPTH    8
#define MAX_NAME_STACK_DEPTH       64
#define MAX_TEXTURE_LEVELS         11
#define MAX_LIGHTS                 16

#define VERTEX_HASH_SIZE 1031

#define MAX_DISPLAY_LISTS 1024
#define OP_BUFFER_MAX_SIZE 512

#define TGL_OFFSET_FILL    0x1
#define TGL_OFFSET_LINE    0x2
#define TGL_OFFSET_POINT   0x4
 
struct GLLight {
  color4f_t ambient;
  color4f_t diffuse;
  color4f_t specular;
  vec4f_t position;	
  vec3f_t spot_direction;
  float spot_exponent;
  float spot_cutoff;
  float attenuation[3];
  /* precomputed values */
  float cos_spot_cutoff;
  vec3f_t norm_spot_direction;
  vec3f_t norm_position;
  /* we use a linked list to know which are the enabled lights */
  int enabled;
  struct GLLight *next,*prev;
};

struct GLMaterial {
  color4f_t emission;
  color4f_t ambient;
  color4f_t diffuse;
  color4f_t specular;
  float shininess;

  /* computed values */
  int shininess_i;
  int do_specular;  
};


struct GLViewport {
  int xmin,ymin,xsize,ysize;
  vec3f_t scale;
  vec3f_t trans;
  int updated;
};

typedef union {
  int op;
  float f;
  int i;
  uint32_t ui;
  void *p;
} GLParam;

  
struct GLParamBuffer {
  GLParam ops[OP_BUFFER_MAX_SIZE];
  struct GLParamBuffer *next;
};

struct GLList {
  GLParamBuffer *first_op_buffer;
  /* TODO: extensions for an hash table or a better allocating scheme */
};

struct GLVertex {
  int edge_flag;
  vec3f_t normal;
  vec4f_t coord;
  vec4f_t tex_coord;
  color4f_t color;
  
  /* computed values */
  vec4f_t ec;                /* eye coordinates */
  vec4f_t pc;                /* coordinates in the normalized volume */
  int clip_code;        /* clip code */
  ZBufferPoint zp;      /* integer coordinates for the rasterization */
};

struct GLImage {
  void *pixmap;
  int xsize,ysize;
};

/* textures */

#define TEXTURE_HASH_TABLE_SIZE 256

struct GLTexture {
  GLImage images[MAX_TEXTURE_LEVELS];
  int handle;
  struct GLTexture *next,*prev;
};


/* shared state */

struct GLSharedState {
  GLList **lists;
  GLTexture **texture_hash_table;
};

struct GLContext;

typedef void (*gl_draw_triangle_func)(struct GLContext *c,
                                      GLVertex *p0,GLVertex *p1,GLVertex *p2);



typedef struct GLSpecBuf {
  int shininess_i;
  int last_used;
  float buf[SPECULAR_BUFFER_SIZE + 1];
  struct GLSpecBuf *next;

  inline void calc_buf(const float shininess) {
    int i;
    float val, inc;
    val = 0.0f;
    inc = 1.0f / SPECULAR_BUFFER_SIZE;
    for (i = 0; i <= SPECULAR_BUFFER_SIZE; i++) {
      buf[i] = pow(val, shininess);
      val += inc;
    }
  }

} GLSpecBuf;
 


/* display context */

struct GLContext {
  /* Z buffer */
  ZBuffer *zb;

  /* lights */
  GLLight lights[MAX_LIGHTS];
  GLLight *first_light;
  color4f_t ambient_light_model;
  int local_light_model;
  int lighting_enabled;
  int light_model_two_side;

  /* materials */
  GLMaterial materials[2];
  int color_material_enabled;
  int current_color_material_mode;
  int current_color_material_type;

  /* textures */
  GLTexture *current_texture;
  int texture_2d_enabled;

  /* shared state */
  GLSharedState shared_state;

  /* current list */
  GLParamBuffer *current_op_buffer;
  int current_op_buffer_index;
  int exec_flag,compile_flag,print_flag;

  /* matrix */

  int matrix_mode;
  mat4_t *matrix_stack[3];
  mat4_t *matrix_stack_ptr[3];
  int matrix_stack_depth_max[3];

  mat4_t matrix_model_view_inv;
  mat4_t matrix_model_projection;
  int matrix_model_projection_updated;
  int matrix_model_projection_no_w_transform; 
  int apply_texture_matrix;

  /* viewport */
  GLViewport viewport;

  /* current state */
  int polygon_mode_back;
  int polygon_mode_front;

  int current_front_face;
  int current_shade_model;
  int current_cull_face;
  int cull_face_enabled;
  int normalize_enabled;
  gl_draw_triangle_func draw_triangle_front,draw_triangle_back;

  /* selection */
  int render_mode;
  uint32_t *select_buffer;
  int select_size;
  uint32_t *select_ptr,*select_hit;
  int select_overflow;
  int select_hits;

  /* names */
  uint32_t name_stack[MAX_NAME_STACK_DEPTH];
  int name_stack_size;

  /* clear */
  float clear_depth;
  color4f_t clear_color;

  /* current vertex state */
  color4f_t current_color;
  uint32_t longcurrent_color[3]; /* precomputed integer color */
  vec4f_t current_normal;
  vec4f_t current_tex_coord;
  int current_edge_flag;

  /* glBegin / glEnd */
  int in_begin;
  int begin_type;
  int vertex_n,vertex_cnt;
  int vertex_max;
  GLVertex *vertex;

  /* opengl 1.1 arrays  */
  float *vertex_array;
  int vertex_array_size;
  int vertex_array_stride;
  float *normal_array;
  int normal_array_stride;
  float *color_array;
  int color_array_size;
  int color_array_stride;
  float *texcoord_array;
  int texcoord_array_size;
  int texcoord_array_stride;
  int client_states;
  
  /* opengl 1.1 polygon offset */
  float offset_factor;
  float offset_units;
  int offset_states;
  
  /* specular buffer. could probably be shared between contexts, 
    but that wouldn't be 100% thread safe */
  GLSpecBuf *specbuf_first;
  int specbuf_used_counter;
  int specbuf_num_buffers;

  /* opaque structure for user's use */
  void *opaque;
  /* resize viewport function */
  int (*gl_resize_viewport)(struct GLContext *c,int *xsize,int *ysize);

  /* depth test */
  int depth_test;

  inline void specbuf_cleanup() {
    /* free all memory used */

    // declared here but never implemented
    // TODO
  } 




  inline GLSpecBuf *specbuf_get_buffer(const int shininess_i,
                                const float shininess) {
    GLSpecBuf *found, *oldest;
    found = oldest = specbuf_first;
    while (found && found->shininess_i != shininess_i) {
      if (found->last_used < oldest->last_used) {
        oldest = found;
      }
      found = found->next;
    }
    if (found) { /* hey, found one! */
      found->last_used = specbuf_used_counter++;
      return found;
    }
    if (oldest == NULL || specbuf_num_buffers < MAX_SPECULAR_BUFFERS) {
      /* create new buffer */
      GLSpecBuf *buf = (GLSpecBuf *)gl_malloc(sizeof(GLSpecBuf));
      if (!buf) gl_fatal_error("could not allocate specular buffer");
      specbuf_num_buffers++;
      buf->next = specbuf_first;
      specbuf_first = buf;
      buf->last_used = specbuf_used_counter++;
      buf->shininess_i = shininess_i;
      buf->calc_buf(shininess);
      return buf;
    }
    /* overwrite the lru buffer */
    /*tgl_trace("overwriting spec buffer :(\n");*/
    oldest->shininess_i = shininess_i;
    oldest->last_used = specbuf_used_counter++;
    oldest->calc_buf(shininess);
    return oldest;
  }

};

extern GLContext *gl_ctx;

void gl_add_op(GLParam *p);

/* clip.c */
void gl_transform_to_viewport(GLContext *c,GLVertex *v);
void gl_draw_triangle(GLContext *c,GLVertex *p0,GLVertex *p1,GLVertex *p2);
void gl_draw_line(GLContext *c,GLVertex *p0,GLVertex *p1);
void gl_draw_point(GLContext *c,GLVertex *p0);

void gl_draw_triangle_point(GLContext *c,
                            GLVertex *p0,GLVertex *p1,GLVertex *p2);
void gl_draw_triangle_line(GLContext *c,
                           GLVertex *p0,GLVertex *p1,GLVertex *p2);
void gl_draw_triangle_fill(GLContext *c,
                           GLVertex *p0,GLVertex *p1,GLVertex *p2);
void gl_draw_triangle_select(GLContext *c,
                             GLVertex *p0,GLVertex *p1,GLVertex *p2);

/* matrix.c */
void gl_print_matrix(const float *m);
/*
void glopLoadIdentity(GLContext *c,GLParam *p);
void glopTranslate(GLContext *c,GLParam *p);*/

/* light.c */
void gl_add_select(GLContext *c,uint32_t zmin,uint32_t zmax);
void gl_enable_disable_light(GLContext *c,int light,int v);
void gl_shade_vertex(GLContext *c,GLVertex *v);

void glInitTextures(GLContext *c);
void glEndTextures(GLContext *c);
GLTexture *alloc_texture(GLContext *c,int h);

/* image_util.c */
void gl_convertRGB_to_5R6G5B(uint16_t *pixmap,uint8_t *rgb,
                             int xsize,int ysize);
void gl_convertRGB_to_8A8R8G8B(uint32_t *pixmap, uint8_t *rgb,
                               int xsize, int ysize);
void gl_resizeImage(uint8_t *dest,int xsize_dest,int ysize_dest,
                    uint8_t *src,int xsize_src,int ysize_src);
void gl_resizeImageNoInterpolate(uint8_t *dest,int xsize_dest,int ysize_dest,
                                 uint8_t *src,int xsize_src,int ysize_src);

GLContext *gl_get_context(void);

void gl_fatal_error(char *format, ...);

 

/* glopXXX functions */

#define ADD_OP(a,b,c) void glop ## a (GLContext *,GLParam *);
#include "opinfo.h"

/* this clip epsilon is needed to avoid some rounding errors after
   several clipping stages */

#define CLIP_EPSILON (1E-5)

static inline int gl_clipcode(float x,float y,float z,float w1)
{
  float w;

  w=w1 * (1.0 + CLIP_EPSILON);
  return (x<-w) |
    ((x>w)<<1) |
    ((y<-w)<<2) |
    ((y>w)<<3) |
    ((z<-w)<<4) | 
    ((z>w)<<5) ;
}

#endif /* _tgl_zgl_h_ */
