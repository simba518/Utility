#include <QMouseEvent>
#include <iostream>
using namespace std;

#include "DragCtrl.h"
using namespace QGLVEXT;

DragCtrl::DragCtrl(pQGLViewerExt_const _viewer):viewer(_viewer){

  begin_drag = false;
  is_dragging = false;
  mouse_button = Qt::LeftButton;
  modify_key = Qt::ControlModifier;
  z_deepth = 0.5f;
}

bool DragCtrl::press (const QMouseEvent *e){

  begin_drag = false;
  is_dragging = false;
  if(viewer != NULL){
	if ((e->button()==mouse_button)&&(e->modifiers() == modify_key)){

	  begin_drag = true;
	  emit startDrag (e->pos().x(),e->pos().y());
	}	
  }
  return begin_drag;
}

bool DragCtrl::release (const QMouseEvent *e){

  const bool perfromed = begin_drag;
  if (begin_drag || is_dragging){
	Vec world_pos = viewer->getWorldCoords(e->pos(),z_deepth);
	emit stopDrag(world_pos[0],world_pos[1],world_pos[2]);
  }
  begin_drag = false;
  is_dragging = false;
  return perfromed;
}

bool DragCtrl::move (const QMouseEvent *e){

  if (begin_drag) {

	initZ_Deepth();
	is_dragging = true;
	begin_drag = false;
	const Vec world_pos = viewer->getWorldCoords(e->pos(),z_deepth);
	emit startDrag(world_pos[0],world_pos[1],world_pos[2]);
  }else if (is_dragging){

	const Vec world_pos = viewer->getWorldCoords(e->pos(),z_deepth);
	emit dragTo(world_pos[0],world_pos[1],world_pos[2]);
  }
  return is_dragging;
}

void DragCtrl::initZ_Deepth(){

  if (drag_hook && viewer){
	double p[3];
	drag_hook->getDragedPoint(p);
	const Vec screen_pos = viewer->getScreenCoords(p[0],p[1],p[2]);
	z_deepth = screen_pos[2];
  }else{
	z_deepth = 0.5f;
  }
}
