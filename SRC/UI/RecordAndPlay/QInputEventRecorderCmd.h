#include <string>
#include <QInputEventRecorder.h>
#include <Log.h>
using namespace std;

// http://algoholic.eu/recording-and-replaying-qt-input-events/

class QInputEventRecorderCmd{

public:
  QInputEventRecorderCmd(QObject *obj=NULL,const bool rec=true,
						 const string f_name="qt_act_recorder.txt"){
	recorder.setObj(obj);
	setCmd(rec,f_name);
  }
  void setCmd(const string cmd,const string f_name){
	setCmd(("replay"==cmd),f_name);
  }
  void setCmd(const bool rec,const string f_name){

	record = rec;
	record_file = QString(f_name.c_str());
	if (!record){
	  recorder.load(f_name.c_str());
	  recorder.replay(1.0);
	}else{
	  record = true;
	  recorder.record();
	}
  }
  ~QInputEventRecorderCmd(){
	if(record)
	  recorder.save(record_file);
  }

private:
  bool record;
  QString record_file;
  QInputEventRecorder recorder;
};
