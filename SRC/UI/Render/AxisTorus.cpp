#include <WindowsHeaders.h>
#include <GL/glew.h>
#include <math.h>
#include "AxisTorus.h"
#include <QGLViewer/config.h>
using namespace QGLVEXT;

void AxisTorus::draw()const{
  
  float scalor[3] = {1.0f ,1.0f ,1.0f};
  if (selected_axis >= 0) {
	scalor[selected_axis%3] = 1.1f;
  }
  drawTorus(0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 0.0f ,scalor[0]);
  drawTorus(1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f ,scalor[1]);
  drawTorus(0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 1.0f ,scalor[2]);
}

void AxisTorus::drawWithNames ()const{

  glFlush();

  glPushName(0);
  drawTorus(0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 0.0f ,1.0f);
  glPopName();

  glPushName(1);
  drawTorus(1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f ,1.0f);
  glPopName();

  glPushName(2);
  drawTorus(0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 1.0f ,1.0f);
  glPopName();
}

void AxisTorus::drawTorus(float x, float y, float z,
						  float r, float g, float b,
						  const float s)const{

  glEnable(GL_COLOR_MATERIAL);
  glColor4f(r,g,b,1.0f);
  glPushMatrix();
  glScalef(s,s,s);
  glRotatef(90.0f, x, y, z);

  const double TWOPI = 2 * M_PI;
  const int numc = 100, numt = 100;

  for (int i = 0; i < numc; i++) {

	glBegin(GL_QUAD_STRIP);
	for (int j = 0; j <= numt; j++) {
	  for (int k = 1; k >= 0; k--) {
		const double s = (i + k) % numc + 0.5;
		const double t = j % numt;
		const double x = (1+0.1 * cos(s * TWOPI/numc))*cos(t*TWOPI/numt);
		const double y = (1+0.1 * cos(s * TWOPI/numc))*sin(t*TWOPI/numt);
		const double z = 0.1 * sin(s * TWOPI / numc);
		glVertex3d(2 * x, 2 * y, 2 * z);
	  }
	}
	glEnd();
  }
  glPopMatrix();
  glDisable(GL_COLOR_MATERIAL);
}
