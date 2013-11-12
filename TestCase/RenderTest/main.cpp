#include <qapplication.h>
#include <QGLViewerExt.h>
#include <MeshRender.h>
#include <eigen3/Eigen/Dense>
using namespace QGLVEXT;
using namespace UTILITY;
using namespace Eigen;

class ObjRender:public SelfRenderEle{
public:
  ObjRender(const Objmesh &obj):_obj(obj){}
  void draw()const{UTILITY::draw(_obj);}
private:
  const Objmesh _obj;
};

class TetRender:public SelfRenderEle{
public:
  TetRender(const TetMesh &tet,const VectorXd &u):_tet(tet),_u(u){}
  TetRender(const TetMesh &tet):_tet(tet){}
  void move(const VectorXd &u){
	assert_eq(u.size(),_tet.nodes().size()*3);
	_u = u;
  }
  void draw()const{
	if (_u.size() > 0){
	  UTILITY::draw(_tet,&_u[0]);
	}else{
	  UTILITY::draw(_tet);
	}
  }
private:
  const TetMesh _tet;
  VectorXd _u; 
};

int main(int argc, char** argv){

  QApplication application(argc,argv);
  QGLViewerExt viewer(NULL);

  const string fn = string(TEST_DATA_DIR)+"beam.obj";
  Objmesh obj;
  obj.load(fn);
  // viewer.addSelfRenderEle(pSelfRenderEle(new ObjRender(obj)));

  const string fn2 = string(TEST_DATA_DIR)+"dino.abq";
  TetMesh tet;
  tet.load(fn2);
  const VectorXd u = VectorXd::Ones(tet.nodes().size()*3)*5.0f;
  viewer.addSelfRenderEle(pSelfRenderEle(new TetRender(tet,u)));
  viewer.addSelfRenderEle(pSelfRenderEle(new TetRender(tet)));

  viewer.setWindowTitle("simpleViewer");
  viewer.show();
  return application.exec();
}
