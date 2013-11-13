#include <string>
#include <QInputEventRecorder.h>
#include <Log.h>
using namespace std;

// http://algoholic.eu/recording-and-replaying-qt-input-events/

#define DEFAULT_RECORD_CMD_FILE "tempt_qt_event_record.b"

class QInputEventRecorderCmd{

public:
  QInputEventRecorderCmd(QObject *obj=NULL,const bool rec=true,
						 const string f_name=DEFAULT_RECORD_CMD_FILE){
	recorder.setObj(obj);
	setCmd(rec,f_name);
  }
  void setCmd(const char *cmd,const string f_name=DEFAULT_RECORD_CMD_FILE){
	INFO_LOG("set command: " << cmd << " " << f_name);
	setCmd((string("record")==string(cmd)),f_name);
  }
  void setCmd(const bool rec,const string f_name=DEFAULT_RECORD_CMD_FILE){

	record = rec;
	record_file = f_name;
	if (!record){
	  INFO_LOG("begin to replay using event file: "<<f_name);
	  recorder.load(f_name.c_str());
	  recorder.replay(1.0);
	}else{
	  INFO_LOG("begin to record using event file: "<<f_name);
	  record = true;
	  recorder.record();
	}
  }
  ~QInputEventRecorderCmd(){
	if(record){
	  INFO_LOG("write events to " << record_file);
	  recorder.save(record_file.c_str());
	}
  }

private:
  bool record;
  string record_file;
  QInputEventRecorder recorder;
};
