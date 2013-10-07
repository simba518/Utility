#include <MeshRender.h>
#include <GL/gl.h>

void UTILITY::draw(const Objmesh& obj){

  const ObjMtl &mtl = obj.getMtl();
  const Eigen::VectorXd &verts = obj.getVerts();
  const Eigen::VectorXi &faces = obj.getFaces();
  const Eigen::VectorXd &norms = obj.getVertNormal();
  assert_eq(verts.size(),norms.size());

  glEnable(GL_SMOOTH);
  glEnable(GL_LIGHTING); 
  glDisable(GL_COLOR_MATERIAL);
  glMaterialfv(GL_FRONT, GL_DIFFUSE, mtl.diffuse);
  glMaterialfv(GL_FRONT, GL_AMBIENT, mtl.ambient);
  glMaterialfv(GL_FRONT, GL_SPECULAR, mtl.specular);
  glMaterialf(GL_FRONT, GL_SHININESS, mtl.shininess);
  glMaterialfv(GL_FRONT, GL_EMISSION, mtl.emission);

  glBegin(GL_TRIANGLES);
  for (int f = 0; f < faces.size(); ++f){
	const int n3 = faces[f]*3;
	assert_in(n3,0,verts.size()-3);
	glNormal3d(norms[n3+0],norms[n3+1],norms[n3+2]);
	glVertex3d(verts[n3+0],verts[n3+1],verts[n3+2]);
  }
  glEnd();
  glDisable(GL_LIGHTING);

}
