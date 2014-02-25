#ifndef _MANIPULATOION_H_
#define _MANIPULATOION_H_

#include <QObject>
#include <QMouseEvent>
#include <QGLViewerExt.h>
#include <eigen3/Eigen/Dense>
#include <assertext.h>
#include <Log.h>
#include <ConTrackBall.h>
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
	LocalframeManipulatoion(pQGLViewerExt v):viewer(v){

	  enabled = false;
	  frame = new ManipulatedFrame();
	  current_pos.setZero();
	  con_track_ball = pConTrackBall(new ConTrackBall(viewer,frame));
	  con_track_ball->setShowConTrackBall(false);
	  connect( frame, SIGNAL(manipulated()), this, SLOT(manipulate()) );
	}
	virtual void setEnable(const bool enable){

	  enabled = enable;
	  if (enabled){

		computeCurrentPosition(current_pos);
		frame->setTranslation(Vec(current_pos[0],current_pos[1],current_pos[2]));
		frame->setOrientation(qglviewer::Quaternion(Vec(1.0,0.0,0.0), 0.0));
		viewer->setManipulatedFrame(frame);

		con_track_ball->translate(current_pos[0],current_pos[1],current_pos[2]);
		const double s = viewer->sceneRadius()/15.0f;
		con_track_ball->scale(s,s,s);
	  }else{
		viewer->setManipulatedFrame(viewer->camera()->frame());
	  }

	  con_track_ball->setEnableConTrackBall(enabled);
	  con_track_ball->setShowConTrackBall(enabled);
	}
	bool isEnable()const{
	  return enabled;
	}
	virtual void draw()const{}
	
	// select
	virtual int totalEleNum()const = 0;
	virtual void drawWithNames()const = 0;
	virtual void select(const vector<int> &sel_ids) = 0;

	// manipulate
	virtual void computeCurrentPosition(Vector3d& pos_xyz)const = 0;
	const Vector3d &currentPosition()const{
	  return current_pos;
	}
	virtual void applyTransform() = 0;

  public slots:
	void manipulate(){
	  computeCurrentPosition(current_pos);
	  con_track_ball->translate(current_pos[0],current_pos[1],current_pos[2]);
	  this->applyTransform();
	}
	
  protected:
	bool enabled;
	pQGLViewerExt viewer;
	ManipulatedFrame* frame;
	pConTrackBall con_track_ball;
	Vector3d current_pos;
  };
  typedef boost::shared_ptr<LocalframeManipulatoion> pLocalframeManipulatoion;

  class LocalframeManipulatoionExt:public LocalframeManipulatoion{

  public:
	LocalframeManipulatoionExt(pQGLViewerExt v):LocalframeManipulatoion(v){}
	virtual void setEnable(const bool enable){

	  LocalframeManipulatoion::setEnable(enable);
	  if (enable){
		computeCurrentPosition(current_pos);
		const Vector3d xc = current_pos;
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
	virtual void computeCurrentPosition(Vector3d &current_pos)const{

		const Matrix<double,3,-1> &xc = this->getCurrentPositions();
		current_pos.setZero();
		for (int i = 0; i < xc.cols(); ++i)
		  current_pos += xc.col(i);
		if (xc.cols() > 0)
		  current_pos *= 1.0f/xc.cols();
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
	  begin_select = false;
  	}
  
  public slots:
  	bool press (QMouseEvent *e){
	  
  	  if ((e->button()==mouse_button)&&(e->modifiers() == modify_key)){
		begin_select = true;
  		manipulate->prepareSelection();
  		viewer->setSelector(manipulate);
  		viewer->select(QRect(e->pos().x(),e->pos().y(),10,10));
  	  }
  	}
  	void select(const vector<int> ids){
	  if (begin_select){
		manipulate->select(ids);
		manipulate->setEnable(ids.size() > 0); 
		begin_select = false;
	  }
  	}

  private:
	bool begin_select;
  	pQGLViewerExt viewer;
  	pLocalframeManipulatoion manipulate;
  	Qt::MouseButton mouse_button;
  	Qt::KeyboardModifiers modify_key;
  };
  typedef boost::shared_ptr<LocalframeManipulatoionCtrl> pLocalframeManipulatoionCtrl;  

}//end of namespace

#endif /*_MANIPULATOION_H_*/
