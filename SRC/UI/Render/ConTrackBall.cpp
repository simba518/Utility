#include <WindowsHeaders.h>
#include <GL/glew.h>
#include <math.h>
#include "ConTrackBall.h"
using namespace qglviewer;
using namespace QGLVEXT;

ConTrackBall::ConTrackBall(pQGLViewerExt viewer):viewer(viewer){

  assert (viewer != NULL);

  p_AxisTorus = pAxisTorus(new AxisTorus());
  constraint = new WorldConstraint();
  grabbed_type = AxisPlaneConstraint::FREE;
  constrained_axi = 0;

  viewer->camera()->frame()->setConstraint(constraint);

  connect( viewer, SIGNAL(mousePressSignal(QMouseEvent *)),this, SLOT(press(QMouseEvent *)) );
  connect( viewer, SIGNAL(mouseReleaseSignal(QMouseEvent *)),this, SLOT(release(QMouseEvent *)) );
  connect( viewer, SIGNAL(selectedIds(vector<int>)),this, SLOT(selectAxises(vector<int>)));
}

void ConTrackBall::selectAxises(const vector<int> sel_group_ids){

  if ( isActive() ){

	if (sel_group_ids.size() <= 0 ){

	  grabbed_type = AxisPlaneConstraint::FORBIDDEN;
	  constrained_axi = -1;
	  p_AxisTorus->selectAxis(constrained_axi);
	}else{

	  grabbed_type = AxisPlaneConstraint::AXIS;
	  constrained_axi = sel_group_ids[0]%3;
	  p_AxisTorus->selectAxis(constrained_axi);
	  viewer->update();
	}
  }
}

void ConTrackBall::checkIfGrabsMouse(int x,int y,const Camera*const camera){

  if ( isFixed() ) {

	viewer->setSelector(p_AxisTorus);
	QRect select_rect(x,y,2,2);
	viewer->select (select_rect);
  }
}

void ConTrackBall::press(QMouseEvent* e){

  if ( isActive() ) {

	constraint->setRotationConstraintType(grabbed_type);
	Vec dir(0.0,0.0,0.0);
	dir[constrained_axi] = 1.0;
	constraint->setRotationConstraintDirection(dir);
  }
}

void ConTrackBall::release(QMouseEvent* e){

  if ( isActive() ){

	grabbed_type = AxisPlaneConstraint::FORBIDDEN;
	constraint->setRotationConstraintType(grabbed_type);
	constrained_axi = -1;
	p_AxisTorus->selectAxis(constrained_axi);
  }
}

void ConTrackBall::setActive(bool ac){

  if (ac){
	constraint->setRotationConstraintType(AxisPlaneConstraint::FORBIDDEN);
	viewer->setMouseTracking(true);
	viewer->addSelfRenderEle(p_AxisTorus);
  }else{
	constraint->setRotationConstraintType(AxisPlaneConstraint::FREE);
	viewer->setMouseTracking(false);
	viewer->removeSelfRenderEle(p_AxisTorus);
  }
}
