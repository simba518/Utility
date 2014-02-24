#ifndef _MANIPULATOION_H_
#define _MANIPULATOION_H_

#include <QObject>
#include <QMouseEvent>
#include <QGLViewerExt.h>
#include <eigen3/Eigen/Dense>
#include <assertext.h>
#include <Log.h>
using namespace qglviewer;
using namespace Eigen;

namespace QGLVEXT{
  
  /**
   * @class LocalframeManipulatoion
   * 
   */
  class LocalframeManipulatoion:public QObject,public SelfRenderEle,public Selectable{

	Q_OBJECT
	
  public:
	LocalframeManipulatoion(QGLViewer *v):viewer(v){

	  frame = new ManipulatedFrame();
	  enabled = false;
	  connect( frame, SIGNAL(manipulated()),this, SLOT(manipulate()) );
	}
	virtual void setEnable(const bool enable){

	  enabled = enable;
	  if (enabled){
		double pos_xyz[3];
		currentPosition(pos_xyz);
		frame->setTranslation(Vec(pos_xyz[0],pos_xyz[1],pos_xyz[2]));
		frame->setOrientation(qglviewer::Quaternion(Vec(1.0,0.0,0.0), 0.0));
		viewer->setManipulatedFrame(frame);
	  }else{
		viewer->setManipulatedFrame(viewer->camera()->frame());
	  }
	}
	bool isEnable()const{
	  return enabled;
	}
	virtual void draw()const{

	  if(enabled){

		double pos_xyz[3];
		currentPosition(pos_xyz);
		glPushMatrix();
		glTranslated(pos_xyz[0],pos_xyz[1],pos_xyz[2]);
		QGLViewer::drawAxis();
		glPopMatrix();
	  }
	}
	
	// select
	virtual int totalEleNum()const = 0;
	virtual void drawWithNames()const = 0;
	virtual void select(const vector<int> &sel_ids) = 0;

	// manipulate
	virtual void currentPosition(double pos_xyz[3])const = 0;
	virtual void applyTransform() = 0;

  public slots:
	void manipulate(){
	  this->applyTransform();
	}
	
  protected:
	QGLViewer *viewer;
	ManipulatedFrame* frame;
	bool enabled;
  };
  typedef boost::shared_ptr<LocalframeManipulatoion> pLocalframeManipulatoion;

  class LocalframeManipulatoionExt:public LocalframeManipulatoion{

  public:
	LocalframeManipulatoionExt(QGLViewer *v):LocalframeManipulatoion(v){}
	virtual void setEnable(const bool enable){

	  LocalframeManipulatoion::setEnable(enable);
	  if (enable){
		double pos_xyz[3];
		currentPosition(pos_xyz);
		const Vector3d xc = Map<Vector3d>(&pos_xyz[0]);
		initial_pos = this->getCurrentPositions();
		for (int i = 0; i < initial_pos.cols(); ++i){
		  initial_pos.col(i) -= xc;
		}
	  }
	}

	// select
	virtual int totalEleNum()const = 0;
	virtual void drawWithNames()const = 0;
	virtual void select(const vector<int> &sel_ids) = 0;
	
	// manipulate
	virtual void currentPosition(double pos_xyz[3])const{

		const Matrix<double,3,-1> &xc = this->getCurrentPositions();
		Vector3d xc_one;
		xc_one.setZero();
		for (int i = 0; i < xc.cols(); ++i)
		  xc_one += xc.col(i);
		if (xc.cols() > 0)
		  xc_one *= 1.0f/xc.cols();
		pos_xyz[0] = xc_one[0];
		pos_xyz[1] = xc_one[1];
		pos_xyz[2] = xc_one[2];
	}
	virtual void applyTransform(){

	  Matrix<double,3,-1> &x = this->getCurrentPositions();
	  assert_eq(x.cols(), initial_pos.cols());
	  assert_eq(x.rows(),3);
	  assert_eq(initial_pos.rows(),3);
	  for (int i = 0; i < initial_pos.cols(); ++i){
		const Vec p0(initial_pos(0,i), initial_pos(1,i), initial_pos(2,i));
		const Vec p = frame->inverseCoordinatesOf(p0);
		x(0,i) = p[0];
		x(1,i) = p[1];
		x(2,i) = p[2];
	  }
	}

	virtual Matrix<double,3,-1> &getCurrentPositions() = 0;
	virtual const Matrix<double,3,-1> &getCurrentPositions()const = 0;

  protected:
	Matrix<double,3,-1> initial_pos;
  };

  class LocalframeManipulatoionCtrl:public QObject{
  
  	Q_OBJECT

  public:
  	LocalframeManipulatoionCtrl(pQGLViewerExt viewer,pLocalframeManipulatoion manipulate):
  	  viewer(viewer), manipulate(manipulate){

  	  mouse_button = Qt::LeftButton;
  	  modify_key = Qt::ControlModifier;
  	  viewer->addSelfRenderEle(manipulate);
  	  connect( viewer, SIGNAL(mousePressSignal(QMouseEvent *)),this, SLOT(press(QMouseEvent *)) );
  	  connect( viewer, SIGNAL(selectedIds(const vector<int> )),this, SLOT(select(const vector<int> )) );
  	}
  
  public slots:
  	bool press (QMouseEvent *e){
  	  if ((e->button()==mouse_button)&&(e->modifiers() == modify_key)){
  		manipulate->prepareSelection();
  		viewer->setSelector(manipulate);
  		viewer->select(QRect(e->pos().x(),e->pos().y(),10,10));
  	  }
  	}
  	void select(const vector<int> ids){
	  manipulate->select(ids);
	  manipulate->setEnable(ids.size() > 0);
  	}

  private:
  	pQGLViewerExt viewer;
  	pLocalframeManipulatoion manipulate;
  	Qt::MouseButton mouse_button;
  	Qt::KeyboardModifiers modify_key;
  };
  typedef boost::shared_ptr<LocalframeManipulatoionCtrl> pLocalframeManipulatoionCtrl;  

}//end of namespace

#endif /*_MANIPULATOION_H_*/
