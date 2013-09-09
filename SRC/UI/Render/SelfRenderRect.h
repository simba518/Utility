#ifndef _SELFRENDERRECT_H_
#define _SELFRENDERRECT_H_

#include <QRect>
#include <QGLViewer/qglviewer.h>
using namespace qglviewer;

#include <SelfRenderEle.h>

namespace QGLVEXT{
  
  class SelfRenderRect: public QRect,public SelfRenderEle{

  public:
	SelfRenderRect(const QGLViewer *p):p_qglviewer(p){
	  setRenderPriority(FOURTH_RENDER);
	}
	void draw()const;
	void clean(){ this->setRect(-1,-1,0,0);}

  private:
	const QGLViewer *p_qglviewer;
  };

  typedef boost::shared_ptr< SelfRenderRect > pSelfRenderRect; 
  
}//end of namespace

#endif /*_SELFRENDERRECT_H_*/