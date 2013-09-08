#include <iostream>
using namespace std;

#include <boost/foreach.hpp>
#include <SelfRenderRect.h>
#include <QMouseEvent>
#include "SelectCtrl.h"

using namespace QGLVEXT;

SelectCtrl::SelectCtrl(pQGLViewerExt viewer,pSelectable selector)
  :viewer(viewer),selector(selector){
 
  assert (viewer != NULL);
  select_rect = pSelfRenderRect( new SelfRenderRect(viewer) );
  select_rect->clean();
  viewer->addSelfRenderEle(select_rect);

  add_mouse_button = Qt::LeftButton;
  add_modify_key = Qt::ShiftModifier;

  rm_mouse_button = Qt::RightButton;
  rm_modify_key = Qt::ShiftModifier;

  sel_status = DO_NOTHING_ON_SEL;
  print_selected_nodes = false;
  enable_op = true;
  begin_select = false;

  createConnections();
}

void SelectCtrl::createConnections(){
  
  // this to viewer
  connect(this ,SIGNAL(rectChanged()), viewer,SLOT(update()));
  connect(this ,SIGNAL(select(QRect)), viewer,SLOT(select(QRect)));

  // viewer to this
  connect( viewer, SIGNAL(mousePressSignal(QMouseEvent *)),this, SLOT(press(QMouseEvent *)) );
  connect( viewer, SIGNAL(mouseReleaseSignal(QMouseEvent *)),this, SLOT(release(QMouseEvent *)) );
  connect( viewer, SIGNAL(mouseMoveSignal(QMouseEvent *)),this, SLOT(move(QMouseEvent *)) );
  connect( viewer, SIGNAL(selectedIds(const vector<int> )),this, SLOT(endSelection(const vector<int> )) );
}

void SelectCtrl::endSelection(const vector<int> sel_ids){

  if (print_selected_nodes) {

	cout <<"select nodes number: "<<sel_ids.size() << endl;
	cout <<"( ";
	BOOST_FOREACH(int ele, sel_ids){
	  cout << ele<< ",";
	}
	cout <<")" << endl;
  }
  if(sel_status == ADD_ELE){
	emit addSelEleMsg(sel_ids);
  }else if(sel_status == REMOVE_ELE){
	emit removeSelEleMsg(sel_ids);
  }
  sel_status = DO_NOTHING_ON_SEL;
}

bool SelectCtrl::press (QMouseEvent *e){

  begin_select = false;
  if (e->modifiers() == add_modify_key&&e->button()==add_mouse_button){

	select_rect->setTopLeft(e->pos());
	select_rect->setBottomRight(e->pos());
	begin_select = true;
	viewer->setSelector(selector);
	sel_status = ADD_ELE;
  }else if (e->modifiers() == rm_modify_key&&e->button()==rm_mouse_button){

	select_rect->setTopLeft(e->pos());
	select_rect->setBottomRight(e->pos());
	begin_select = true;
	viewer->setSelector(selector);
	sel_status = REMOVE_ELE;
  }

  return begin_select;
}

bool SelectCtrl::release (QMouseEvent *e){

  bool perferm = false;
  if(begin_select){

	select_rect->normalized();
	emit select(*select_rect);
	select_rect->clean();
	begin_select = false;
	emit rectChanged();
  }
  return perferm;
}

bool SelectCtrl::move (QMouseEvent *e){

  bool perferm = false;
  if(begin_select){

	select_rect->setBottomRight(e->pos());
	emit rectChanged();
  }
  return perferm;
}
