#include "msghandling.hpp"
#include "zgl.hpp"

void glopMaterial(GLContext *c, GLParam *p) {
  int mode = p[1].i;
  int type = p[2].i;
  int i;
  GLMaterial *m;

  if (mode == GL_FRONT_AND_BACK) {
    p[1].i = GL_FRONT;
    glopMaterial(c, p);
    mode = GL_BACK;
  }
  if (mode == GL_FRONT)
    m = &c->materials[0];
  else
    m = &c->materials[1];

  switch (type) {
    case GL_EMISSION:
      m->emission = rgba_t(p[3 + 0].f, p[3 + 1].f, p[3 + 2].f, p[3 + 3].f);
      break;
    case GL_AMBIENT:
      m->ambient = rgba_t(p[3 + 0].f, p[3 + 1].f, p[3 + 2].f, p[3 + 3].f);
      break;
    case GL_DIFFUSE:
      m->diffuse = rgba_t(p[3 + 0].f, p[3 + 1].f, p[3 + 2].f, p[3 + 3].f);
      break;
    case GL_SPECULAR:
      m->specular = rgba_t(p[3 + 0].f, p[3 + 1].f, p[3 + 2].f, p[3 + 3].f);
      break;
    case GL_SHININESS:
      m->shininess = p[3].f;
      m->shininess_i = (m->shininess / 128.0f) * SPECULAR_BUFFER_RESOLUTION;
      break;
    case GL_AMBIENT_AND_DIFFUSE:
      m->diffuse = rgba_t(p[3 + 0].f, p[3 + 1].f, p[3 + 2].f, p[3 + 3].f);
      m->ambient = rgba_t(p[3 + 0].f, p[3 + 1].f, p[3 + 2].f, p[3 + 3].f);
      break;
    default:
      assert(0);
  }
}

void glopColorMaterial(GLContext *c, GLParam *p) {
  int mode = p[1].i;
  int type = p[2].i;

  c->current_color_material_mode = mode;
  c->current_color_material_type = type;
}

void glopLight(GLContext *c, GLParam *p) {
  int light = p[1].i;
  int type = p[2].i;
  float v[4];
  GLLight *l;
  int i;

  assert(light >= GL_LIGHT0 && light < GL_LIGHT0 + MAX_LIGHTS);

  l = &c->lights[light - GL_LIGHT0];
  for (i = 0; i < 4; i++) v[i] = p[3 + i].f;

  switch (type) {
    case GL_AMBIENT:
      l->ambient = rgba_t(v[0], v[1], v[2], v[3]);
      break;
    case GL_DIFFUSE:
      l->diffuse = rgba_t(v[0], v[1], v[2], v[3]);
      break;
    case GL_SPECULAR:
      l->specular = rgba_t(v[0], v[1], v[2], v[3]);
      break;
    case GL_POSITION: {
      vec4_t pos;
      vec4_t _v = vec4_t(v[0], v[1], v[2], v[3]);
      pos = c->matrix_stack_ptr[0]->Mulvec4_t(_v);

      l->position = pos;

      if (l->position.w == 0) {
        l->norm_position.x = pos.x;
        l->norm_position.y = pos.y;
        l->norm_position.z = pos.z;

        l->norm_position.Norm();
      }
    } break;
    case GL_SPOT_DIRECTION:
      l->spot_direction = vec3_t(v[0], v[1], v[2]);
      l->norm_spot_direction = vec3_t(v[0], v[1], v[2]);
      l->norm_spot_direction.Norm();
      break;
    case GL_SPOT_EXPONENT:
      l->spot_exponent = v[0];
      break;
    case GL_SPOT_CUTOFF: {
      float a = v[0];
      assert(a == 180 || (a >= 0 && a <= 90));
      l->spot_cutoff = a;
      if (a != 180) l->cos_spot_cutoff = cos(a * M_PI / 180.0);
    } break;
    case GL_CONSTANT_ATTENUATION:
      l->attenuation[0] = v[0];
      break;
    case GL_LINEAR_ATTENUATION:
      l->attenuation[1] = v[0];
      break;
    case GL_QUADRATIC_ATTENUATION:
      l->attenuation[2] = v[0];
      break;
    default:
      assert(0);
  }
}

void glopLightModel(GLContext *c, GLParam *p) {
  int pname = p[1].i;
  int i;

  switch (pname) {
    case GL_LIGHT_MODEL_AMBIENT:
      c->ambient_light_model =
          rgba_t(p[2 + 0].f, p[2 + 1].f, p[2 + 2].f, p[2 + 3].f);
      break;
    case GL_LIGHT_MODEL_LOCAL_VIEWER:
      c->local_light_model = (int)p[2].f;
      break;
    case GL_LIGHT_MODEL_TWO_SIDE:
      c->light_model_two_side = (int)p[2].f;
      break;
    default:
      tgl_warning("glopLightModel: illegal pname: 0x%x\n", pname);
      // assert(0);
      break;
  }
}

static inline float clampf(float a, float min, float max) {
  if (a < min)
    return min;
  else if (a > max)
    return max;
  else
    return a;
}

void gl_enable_disable_light(GLContext *c, int light, int v) {
  GLLight *l = &c->lights[light];
  if (v && !l->enabled) {
    l->enabled = 1;
    l->next = c->first_light;
    c->first_light = l;
    l->prev = NULL;
  } else if (!v && l->enabled) {
    l->enabled = 0;
    if (l->prev == NULL)
      c->first_light = l->next;
    else
      l->prev->next = l->next;
    if (l->next != NULL) l->next->prev = l->prev;
  }
}

/* non optimized lightening model */
void gl_shade_vertex(GLContext *c, GLVertex *v) {
  float R, G, B, A;
  GLMaterial *m;
  GLLight *l;
  vec3_t n, s, d;
  float dist, tmp, att, dot, dot_spot, dot_spec;
  int twoside = c->light_model_two_side;

  m = &c->materials[0];

  n.x = v->normal.x;
  n.y = v->normal.y;
  n.z = v->normal.z;

  R = m->emission.R + m->ambient.R * c->ambient_light_model.R;
  G = m->emission.G + m->ambient.G * c->ambient_light_model.G;
  B = m->emission.B + m->ambient.B * c->ambient_light_model.B;
  A = clampf(m->diffuse.A, 0, 1);

  for (l = c->first_light; l != NULL; l = l->next) {
    float lR, lB, lG;

    /* ambient */
    lR = l->ambient.R * m->ambient.R;
    lG = l->ambient.G * m->ambient.G;
    lB = l->ambient.B * m->ambient.B;

    if (l->position.w == 0) {
      /* light at infinity */
      d.x = l->norm_position.x;
      d.y = l->norm_position.y;
      d.z = l->norm_position.z;
      att = 1;
    } else {
      /* distance attenuation */
      d.x = l->position.x - v->ec.x;
      d.y = l->position.y - v->ec.y;
      d.z = l->position.z - v->ec.z;
      dist = sqrt(d.x * d.x + d.y * d.y + d.z * d.z);
      if (dist > 1E-10f) {
        tmp = 1 / dist;
        d.x *= tmp;
        d.y *= tmp;
        d.z *= tmp;
      }
      att = 1.0f / (l->attenuation[0] +
                    dist * (l->attenuation[1] + dist * l->attenuation[2]));
    }
    dot = d.x * n.x + d.y * n.y + d.z * n.z;
    if (twoside && dot < 0) dot = -dot;
    if (dot > 0) {
      /* diffuse light */
      lR += dot * l->diffuse.R * m->diffuse.R;
      lG += dot * l->diffuse.G * m->diffuse.G;
      lB += dot * l->diffuse.B * m->diffuse.B;

      /* spot light */
      if (l->spot_cutoff != 180) {
        dot_spot =
            -(d.x * l->norm_spot_direction.x + d.y * l->norm_spot_direction.y +
              d.z * l->norm_spot_direction.z);
        if (twoside && dot_spot < 0) dot_spot = -dot_spot;
        if (dot_spot < l->cos_spot_cutoff) {
          /* no contribution */
          continue;
        } else {
          /* TODO: optimize */
          if (l->spot_exponent > 0) {
            att = att * pow(dot_spot, l->spot_exponent);
          }
        }
      }

      /* specular light */

      if (c->local_light_model) {
        vec3_t vcoord;
        vcoord.x = v->ec.x;
        vcoord.y = v->ec.y;
        vcoord.z = v->ec.z;
        vcoord.Norm();
        s.x = d.x - vcoord.x;
        s.y = d.y - vcoord.x;
        s.z = d.z - vcoord.x;
      } else {
        s.x = d.x;
        s.y = d.y;
        s.z = d.z + 1.0;
      }
      dot_spec = n.x * s.x + n.y * s.y + n.z * s.z;
      if (twoside && dot_spec < 0) dot_spec = -dot_spec;
      if (dot_spec > 0) {
        GLSpecBuf *specbuf;
        int idx;
        tmp = sqrt(s.x * s.x + s.y * s.y + s.z * s.z);
        if (tmp > 1E-3) {
          dot_spec = dot_spec / tmp;
        }

        /* TODO: optimize */
        /* testing specular buffer code */
        /* dot_spec= pow(dot_spec,m->shininess);*/
        specbuf = c->specbuf_get_buffer(m->shininess_i, m->shininess);
        idx = (int)(dot_spec * SPECULAR_BUFFER_SIZE);
        if (idx > SPECULAR_BUFFER_SIZE) idx = SPECULAR_BUFFER_SIZE;
        dot_spec = specbuf->buf[idx];
        lR += dot_spec * l->specular.R * m->specular.R;
        lG += dot_spec * l->specular.G * m->specular.G;
        lB += dot_spec * l->specular.B * m->specular.B;
      }
    }

    R += att * lR;
    G += att * lG;
    B += att * lB;
  }

  v->color.R = clampf(R, 0, 1);
  v->color.G = clampf(G, 0, 1);
  v->color.B = clampf(B, 0, 1);
  v->color.A = A;
}
