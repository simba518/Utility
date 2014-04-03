#include <string>
#include <iostream>
#include <sstream>
#include "QInputEventRecorderWidget.h"
using namespace std;
using namespace UTILITY;

#define TOSTR(input) boost::lexical_cast<std::string>(input)

QInputEventRecorderWidget::QInputEventRecorderWidget(QWidget *parent,Qt::WFlags flags)
:QDialog(parent,flags){

  dialog.setupUi(this);
  createConnections();
  updateTimeOffset();
}

void QInputEventRecorderWidget::createConnections(){

  connect(dialog.pushButton_load, SIGNAL(clicked()),this,SLOT(load()));
  connect(dialog.pushButton_save, SIGNAL(clicked()),this,SLOT(save()));
  connect(dialog.pushButton_cut,SIGNAL(clicked()),&recorder,SLOT(removeTheUnplayedEvents()));
  connect(dialog.pushButton_clear,SIGNAL(clicked()),&recorder,SLOT(clear()));

  connect(dialog.radioButton_pause, SIGNAL(clicked()), &recorder, SLOT(stop()));
  connect(dialog.radioButton_play, SIGNAL(clicked()), &recorder, SLOT(replay()));
  connect(dialog.radioButton_record, SIGNAL(clicked()), &recorder, SLOT(record()));

  connect(&recorder,SIGNAL(isPaused(bool)),dialog.radioButton_pause,SLOT(setChecked(bool)));
  connect(&recorder,SIGNAL(updateTimeOffset()),this,SLOT(updateTimeOffset()));
}

void QInputEventRecorderWidget::updateTimeOffset(){
  
  const int total = recorder.totalTimeOffset();
  const int current = recorder.currentTimeOffset();

  dialog.progressBar->setRange(0,total>0?total:1);
  dialog.progressBar->setValue(current);

  stringstream ss;
  ss.precision(5);
  ss << (current/1000.0f) << " / " << (total/1000.0f) << "s";
  const string ssd = ss.str(); 
  dialog.current_total->setText(ssd.c_str());
}
