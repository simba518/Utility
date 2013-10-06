#include <iostream>
#include <map>
#include <string>
using namespace std;

typedef struct{
  int num;
  string name;
}Student;

int main(int argc, char *argv[]){

  map<int,string> stu;
  stu[1] = "a;dkfjasd;lkfj";
  map<int,string>::const_iterator it;
  it = stu.find(1);
  if (it != stu.end()){
	cout << it->first << endl;
	cout << it->second << endl;
  }

  // Student s[] = {
  // 	{100,"hh"},
  // 	{101,"hh"},
  // 	{102,"haa"}
  // };
  
  return 0;
}
