#ifndef _DRAG3D_H_
#define _DRAG3D_H_

#include <QObject>
#include <QGLViewerExt.h>
using namespace qglviewer;

namespace QGLVEXT{

  /**
   * @class DragHook return the position of the point that is dragged, which is
   * used to set the z-deepth for the drag operation.
   * 
   */
  class DragHook{
	
  public:
	virtual void getDragedPoint(const double point[3])const = 0;

  };
  typedef boost::shared_ptr<DragHook> pDragHook;

  /**
   * @class DragCtrl control the drag operation.
   */
  class DragCtrl:public QObject{

	Q_OBJECT
	
  public:
	DragCtrl(pQGLViewerExt_const viewer);
	void setKeyMouse(Qt::KeyboardModifiers key,Qt::MouseButton m_button){
	  this->modify_key = key;
	  this->mouse_button = m_button;
	}
	bool beginDrag()const{
	  return begin_drag;
	}
	void setDragHook(pDragHook drag_hook){
	  this->drag_hook = drag_hook;
	}
	void initZ_Deepth();

  protected:
	bool press (const QMouseEvent *e);
	bool release (const QMouseEvent *e);
	bool move (const QMouseEvent *e);
	
  signals:
	void startDrag (int screen_x,int screen_y);
	void startDrag (double x,double y,double z);
	void dragTo (double x,double y,double z);
	void stopDrag (double x,double y,double z);
	
  private:
	pQGLViewerExt_const viewer;
	pDragHook drag_hook;

	/**
	 * the mouse button and modify key that trags the draging operation, could
	 * be set by calling setKeyMouse(), and the default value is left+Ctrl.
	 */
	Qt::MouseButton mouse_button;
	Qt::KeyboardModifiers modify_key;

	double z_deepth;
	bool begin_drag;
	bool is_dragging;
  };

  typedef boost::shared_ptr< DragCtrl > pDragCtrl;
  
}//end of namespace

#endif /*_DRAG3D_H_*/
