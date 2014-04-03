#ifndef _QINPUTEVENTRECORDERWIDGET_H_
#define _QINPUTEVENTRECORDERWIDGET_H_

#include <QFileDialog>
#include <QInputEventRecorder.h>
#include "ui_record_replay.h"

namespace UTILITY{
  
  /**
   * @class QInputEventRecorderWidget
   * 
   */
  class QInputEventRecorderWidget:public QDialog{

	Q_OBJECT
	
  public:
	QInputEventRecorderWidget(QWidget *parent=0,Qt::WFlags flags=0);
	void setObj(QObject* widget){
	  recorder.setObj(widget);
	}

  protected:
	void createConnections();
	
  protected slots:
	void load(){
	  recorder.stop();
	  QString filename = QFileDialog::getOpenFileName(NULL);
	  if (!filename.isEmpty()){
		recorder.load(filename);
	  }
	}
	void save(){
	  recorder.stop();
	  const QString filename = QFileDialog::getSaveFileName(NULL);
	  if (!filename.isEmpty()){
		recorder.save(filename);
	  }
	}
	void updateTimeOffset();

  private:
	QInputEventRecorder recorder;
	Ui_record_replay_dialog dialog;
  };
  
}//end of namespace

#endif /*_QINPUTEVENTRECORDERWIDGET_H_*/
