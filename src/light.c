#include "zgl.hpp"
#include "msghandling.h"

void glopMaterial(GLContext *c,GLParam *p)
{
  int mode=p[1].i;
  int type=p[2].i;
  int i;
  GLMaterial *m;

  if (mode == GL_FRONT_AND_BACK) {
    p[1].i=GL_FRONT;
    glopMaterial(c,p);
    mode=GL_BACK;
  }
  if (mode == GL_FRONT) m=&c->materials[0];
  else m=&c->materials[1];
  
  switch(type) {
  case GL_EMISSION: 
      m->emission=COLOR4(p[3+0].f,p[3+1].f,p[3+2].f,p[3+3].f);
    break;
  case GL_AMBIENT:
      m->ambient=COLOR4(p[3+0].f,p[3+1].f,p[3+2].f,p[3+3].f);
    break;
  case GL_DIFFUSE:
      m->diffuse=COLOR4(p[3+0].f,p[3+1].f,p[3+2].f,p[3+3].f);
    break;
  case GL_SPECULAR:
      m->specular=COLOR4(p[3+0].f,p[3+1].f,p[3+2].f,p[3+3].f);
    break;
  case GL_SHININESS:
    m->shininess=p[3].f;
    m->shininess_i = (m->shininess/128.0f)*SPECULAR_BUFFER_RESOLUTION;
    break;
  case GL_AMBIENT_AND_DIFFUSE:  
    m->diffuse=COLOR4(p[3+0].f,p[3+1].f,p[3+2].f,p[3+3].f); 
    m->ambient=COLOR4(p[3+0].f,p[3+1].f,p[3+2].f,p[3+3].f);
    break;
  default:
    assert(0);
  }
}

void glopColorMaterial(GLContext *c,GLParam *p)
{
  int mode=p[1].i;
  int type=p[2].i;

  c->current_color_material_mode=mode;
  c->current_color_material_type=type;
}

void glopLight(GLContext *c,GLParam *p)
{
  int light=p[1].i;
  int type=p[2].i;
  float v[4];
  GLLight *l;
  int i;
  
  assert(light >= GL_LIGHT0 && light < GL_LIGHT0+MAX_LIGHTS );

  l=&c->lights[light-GL_LIGHT0];
  for(i=0;i<4;i++) v[i]=p[3+i].f;

  switch(type) {
  case GL_AMBIENT:
    l->ambient=COLOR4(v[0],v[1],v[2],v[3]);
    break;
  case GL_DIFFUSE:
    l->diffuse=COLOR4(v[0],v[1],v[2],v[3]);
    break;
  case GL_SPECULAR:
    l->specular=COLOR4(v[0],v[1],v[2],v[3]);
    break;
  case GL_POSITION:
    {
      V4 pos;
      V4 _v=V4(v[0],v[1],v[2],v[3]); 
      pos =c->matrix_stack_ptr[0]->MulV4(_v);

      l->position=pos;

      if (l->position.W == 0) {
        l->norm_position.X=pos.X;
        l->norm_position.Y=pos.Y;
        l->norm_position.Z=pos.Z;
        
        l->norm_position.Norm();
      }
    }
    break;
  case GL_SPOT_DIRECTION:  
      l->spot_direction=V3(v[0],v[1],v[2]); 
      l->norm_spot_direction=V3(v[0],v[1],v[2]); 
    l->norm_spot_direction.Norm();
    break;
  case GL_SPOT_EXPONENT:
    l->spot_exponent=v[0];
    break;
  case GL_SPOT_CUTOFF:
    {
      float a=v[0];
      assert(a == 180 || (a>=0 && a<=90));
      l->spot_cutoff=a;
      if (a != 180) l->cos_spot_cutoff=cos(a * M_PI / 180.0);
    }
    break;
  case GL_CONSTANT_ATTENUATION:
    l->attenuation[0]=v[0];
    break;
  case GL_LINEAR_ATTENUATION:
    l->attenuation[1]=v[0];
    break;
  case GL_QUADRATIC_ATTENUATION:
    l->attenuation[2]=v[0];
    break;
  default:
    assert(0);
  }
}
  

void glopLightModel(GLContext *c,GLParam *p)
{
  int pname=p[1].i;
  int i;

  switch(pname) {
  case GL_LIGHT_MODEL_AMBIENT: 
      c->ambient_light_model=COLOR4(p[2 + 0].f,p[2 + 1].f,p[2 + 2].f,p[2 + 3].f);
    break;
  case GL_LIGHT_MODEL_LOCAL_VIEWER:
    c->local_light_model=(int)p[2].f;
    break;
  case GL_LIGHT_MODEL_TWO_SIDE:
    c->light_model_two_side = (int)p[2].f;
    break;
  default:
    tgl_warning("glopLightModel: illegal pname: 0x%x\n", pname);
    //assert(0);
    break;
  }
}


static inline float clampf(float a,float min,float max)
{
  if (a<min) return min;
  else if (a>max) return max;
  else return a;
}

void gl_enable_disable_light(GLContext *c,int light,int v)
{
  GLLight *l=&c->lights[light];
  if (v && !l->enabled) {
    l->enabled=1;
    l->next=c->first_light;
    c->first_light=l;
    l->prev=NULL;
  } else if (!v && l->enabled) {
    l->enabled=0;
    if (l->prev == NULL) c->first_light=l->next;
    else l->prev->next=l->next;
    if (l->next != NULL) l->next->prev=l->prev;
  }
}

/* non optimized lightening model */
void gl_shade_vertex(GLContext *c,GLVertex *v)
{
  float R,G,B,A;
  GLMaterial *m;
  GLLight *l;
  V3 n,s,d;
  float dist,tmp,att,dot,dot_spot,dot_spec;
  int twoside = c->light_model_two_side;

  m=&c->materials[0];

  n.X=v->normal.X;
  n.Y=v->normal.Y;
  n.Z=v->normal.Z;

  R=m->emission.R+m->ambient.R*c->ambient_light_model.R;
  G=m->emission.G+m->ambient.G*c->ambient_light_model.G;
  B=m->emission.B+m->ambient.B*c->ambient_light_model.B;
  A=clampf(m->diffuse.A,0,1);

  for(l=c->first_light;l!=NULL;l=l->next) {
    float lR,lB,lG;
    
    /* ambient */
    lR=l->ambient.R * m->ambient.R;
    lG=l->ambient.G * m->ambient.G;
    lB=l->ambient.B * m->ambient.B;

    if (l->position.W == 0) {
      /* light at infinity */
      d.X=l->norm_position.X;
      d.Y=l->norm_position.Y;
      d.Z=l->norm_position.Z;
      att=1;
    } else {
      /* distance attenuation */
      d.X=l->position.X-v->ec.X;
      d.Y=l->position.Y-v->ec.Y;
      d.Z=l->position.Z-v->ec.Z;
      dist=sqrt(d.X*d.X+d.Y*d.Y+d.Z*d.Z);
      if (dist>1E-10f) {
        tmp=1/dist;
        d.X*=tmp;
        d.Y*=tmp;
        d.Z*=tmp;
      }
      att=1.0f/(l->attenuation[0]+dist*(l->attenuation[1]+
				     dist*l->attenuation[2]));
    }
    dot=d.X*n.X+d.Y*n.Y+d.Z*n.Z;
    if (twoside && dot < 0) dot = -dot;
    if (dot>0) {
      /* diffuse light */
      lR+=dot * l->diffuse.R * m->diffuse.R;
      lG+=dot * l->diffuse.G * m->diffuse.G;
      lB+=dot * l->diffuse.B * m->diffuse.B;

      /* spot light */
      if (l->spot_cutoff != 180) {
        dot_spot=-(d.X*l->norm_spot_direction.X+
                   d.Y*l->norm_spot_direction.Y+
                   d.Z*l->norm_spot_direction.Z);
        if (twoside && dot_spot < 0) dot_spot = -dot_spot;
        if (dot_spot < l->cos_spot_cutoff) {
          /* no contribution */
          continue;
        } else {
          /* TODO: optimize */
          if (l->spot_exponent > 0) {
            att=att*pow(dot_spot,l->spot_exponent);
          }
        }
      }

      /* specular light */
      
      if (c->local_light_model) {
        V3 vcoord;
        vcoord.X=v->ec.X;
        vcoord.Y=v->ec.Y;
        vcoord.Z=v->ec.Z;
        vcoord.Norm();
        s.X=d.X-vcoord.X;
        s.Y=d.Y-vcoord.X;
        s.Z=d.Z-vcoord.X;
      } else {
        s.X=d.X;
        s.Y=d.Y;
        s.Z=d.Z+1.0;
      }
      dot_spec=n.X*s.X+n.Y*s.Y+n.Z*s.Z;
      if (twoside && dot_spec < 0) dot_spec = -dot_spec;
      if (dot_spec>0) {
        GLSpecBuf *specbuf;
        int idx;
        tmp=sqrt(s.X*s.X+s.Y*s.Y+s.Z*s.Z);
        if (tmp > 1E-3) {
          dot_spec=dot_spec / tmp;
        }
      
        /* TODO: optimize */
        /* testing specular buffer code */
        /* dot_spec= pow(dot_spec,m->shininess);*/
        specbuf = specbuf_get_buffer(c, m->shininess_i, m->shininess);
        idx = (int)(dot_spec*SPECULAR_BUFFER_SIZE);
        if (idx > SPECULAR_BUFFER_SIZE) idx = SPECULAR_BUFFER_SIZE;
        dot_spec = specbuf->buf[idx];
        lR+=dot_spec * l->specular.R * m->specular.R;
        lG+=dot_spec * l->specular.G * m->specular.G;
        lB+=dot_spec * l->specular.B * m->specular.B;
      }
    }

  
    R+=att * lR;
    G+=att * lG;
    B+=att * lB;
  }

  v->color.R=clampf(R,0,1);
  v->color.G=clampf(G,0,1);
  v->color.B=clampf(B,0,1);
  v->color.A=A;
}

