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

int main(int argc, char** argv){

  QApplication application(argc,argv);
  QGLViewerExt viewer(NULL);

  // const string fn="/home/simba/Workspace/AnimationEditor/Data/bunny/bunny.obj";
  // const string fn="/home/simba/Workspace/AnimationEditor/Data/beam/beam.obj";
  const string fn = "./TestCase/TestData/beam.obj";
  Objmesh obj;
  obj.load(fn);
  viewer.addSelfRenderEle(pSelfRenderEle(new ObjRender(obj)));

  viewer.setWindowTitle("simpleViewer");
  // viewer.toggleDrawLights();
  viewer.show();
  return application.exec();
}
