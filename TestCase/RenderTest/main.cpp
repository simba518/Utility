#include <qapplication.h>
#include <QGLViewerExt.h>
#include <QInputEventRecorderCmd.h>
#include "QtEventRecorder.h"
#include "SimpleMeshRender.h"
using namespace QGLVEXT;
using namespace UTILITY;
using namespace Eigen;

int main(int argc, char** argv){

  QApplication application(argc,argv);
  QGLViewerExt viewer(NULL);
  QInputEventRecorderCmd recorder(NULL,false);

  {
	// const string fn = string(TEST_DATA_DIR)+"beam.obj";
	// Objmesh obj;
	// obj.load(fn);
	// viewer.addSelfRenderEle(pSelfRenderEle(new ObjRender(obj)));

	const string fn2 = string(TEST_DATA_DIR)+"dino.abq";
	TetMesh tet;
	tet.load(fn2);
	const VectorXd u = VectorXd::Ones(tet.nodes().size()*3)*5.0f;
	viewer.addSelfRenderEle(pSelfRenderEle(new TetRender(tet,u)));
	viewer.addSelfRenderEle(pSelfRenderEle(new TetRender(tet)));
  }

  viewer.setWindowTitle("simpleViewer");
  viewer.show();
  return application.exec();
}
