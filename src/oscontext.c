#include <GL/gl.h>
#include <GL/oscontext.h>
#include <assert.h>
#include <stdlib.h>

#include "zbuffer.hpp"
#include "zgl.hpp"

static int buffercnt = 0;

ostgl_context *ostgl_create_context(const int xsize, const int ysize,
                                    const int depth, void **framebuffers,
                                    const int numbuffers) {
  ostgl_context *context;
  int i;
  ZBuffer *zb;

  assert(depth == 16); /* support for other depths must include bpp
                          convertion */
  assert(numbuffers >= 1);

  context = (ostgl_context *)gl_malloc(sizeof(ostgl_context));
  assert(context);
  context->zbs = (void **)gl_malloc(sizeof(void *) * numbuffers);
  context->framebuffers = (void **)gl_malloc(sizeof(void *) * numbuffers);

  assert(context->zbs != NULL && context->framebuffers != NULL);

  for (i = 0; i < numbuffers; i++) {
    context->framebuffers[i] = framebuffers[i];
    zb = ZBuffer::open(xsize, ysize, 0, NULL, NULL, framebuffers[i]);
    if (zb == NULL) {
      fprintf(stderr, "Error while initializing Z buffer\n");
      exit(1);
    }
    context->zbs[i] = zb;
  }
  if (++buffercnt == 1) {
    glInit(context->zbs[0]);
  }
  context->xsize = xsize;
  context->ysize = ysize;
  context->numbuffers = numbuffers;
  return context;
}

void ostgl_delete_context(ostgl_context *context) {
  int i;
  for (i = 0; i < context->numbuffers; i++) {
    ZBuffer::close((ZBuffer *)context->zbs[i]);
  }
  gl_free(context->zbs);
  gl_free(context->framebuffers);
  gl_free(context);

  if (--buffercnt == 0) {
    glClose();
  }
}

void ostgl_make_current(ostgl_context *oscontext, const int idx) {
  GLContext *context = gl_get_context();
  assert(idx < oscontext->numbuffers);
  context->zb = (ZBuffer *)(oscontext->zbs[idx]);
}

void ostgl_resize(ostgl_context *context, const int xsize, const int ysize,
                  void **framebuffers) {
  int i;
  for (i = 0; i < context->numbuffers; i++) {
    ((ZBuffer *)context->zbs[i])->resize(framebuffers[i], xsize, ysize);
  }
}
