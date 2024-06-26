#include "zgl.hpp"

void glopClearColor(GLContext *c, GLParam *p) {
  c->clear_color.r = p[1].f;
  c->clear_color.g = p[2].f;
  c->clear_color.b = p[3].f;
  c->clear_color.a = p[4].f;
}
void glopClearDepth(GLContext *c, GLParam *p) { c->clear_depth = p[1].f; }

void glopClear(GLContext *c, GLParam *p) {
  int mask = p[1].i;
  int z = 0;
  int r = (int)(c->clear_color.r * 65535);
  int g = (int)(c->clear_color.g * 65535);
  int b = (int)(c->clear_color.b * 65535);

  /* TODO : correct value of Z */

  c->zb->clear(mask & GL_DEPTH_BUFFER_BIT, z, mask & GL_COLOR_BUFFER_BIT, r, g,
               b);
}
