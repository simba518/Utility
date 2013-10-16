#include <MeshRender.h>
#include <GL/gl.h>

void UTILITY::draw(const Objmesh& obj){

  const ObjMtl &mtl = obj.getMtl();
  const Eigen::VectorXd &verts = obj.getVerts();
  const Eigen::VectorXd &norms = obj.getVertNormal();
  const Eigen::VectorXi &faces = obj.getFaces();
  const Eigen::VectorXi &normIndex = obj.getNormalIndex();

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
	const int v3 = faces[f]*3;
	const int n3 = normIndex[f]*3;
	assert_in(v3,0,verts.size()-3);
	assert_in(n3,0,norms.size()-3);
	glNormal3d(norms[n3+0],norms[n3+1],norms[n3+2]);
	glVertex3d(verts[v3+0],verts[v3+1],verts[v3+2]);
  }
  glEnd();
  glDisable(GL_LIGHTING);

}
