#ifndef _CONTRACKBALL_H_
#define _CONTRACKBALL_H_

#include <boost/shared_ptr.hpp>
#include <QGLViewerExt.h>
#include <SelfRenderEle.h>
#include <AxisTorus.h>

namespace QGLVEXT{

  /**
   * @class ConTrackBall constrained track ball. the rotation will be
   * constrained to selected axis.
   * 
   */
  class ConTrackBall: public QObject, public MouseGrabber{

	Q_OBJECT

  public:
  	ConTrackBall(pQGLViewerExt viewer);
  	void checkIfGrabsMouse(int x, int y, const Camera* const camera);

  	void setActive(bool ac);
  	bool isActive()const{
	  return (constraint->rotationConstraintType() != AxisPlaneConstraint::FREE);
	}
  	bool isFixed()const{
	  return (AxisPlaneConstraint::FORBIDDEN==constraint->rotationConstraintType());
	}

  public slots:
  	void selectAxises(const vector<int> sel_group_ids);
  	void press(QMouseEvent* e);
  	void release(QMouseEvent* e);
  	void toggleActive(){
	  setActive( !isActive() );
  	}

  private:
  	pQGLViewerExt viewer;
  	pAxisTorus p_AxisTorus;
  	AxisPlaneConstraint* constraint;
  	AxisPlaneConstraint::Type grabbed_type;
  	short constrained_axi; // 0: x, 1: y, 2: z, -1:none.
  };
  
  typedef boost::shared_ptr<ConTrackBall> pConTrackBall;
  
}//end of namespace

#endif /*_CONTRACKBALL_H_*/
