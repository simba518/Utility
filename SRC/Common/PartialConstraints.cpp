#include <AuxTools.h>
#include <map>
#include <boost/foreach.hpp>
#include <assertext.h>
#include "PartialConstraints.h"
using namespace std;
using namespace UTILITY;

bool PartialConstraints::write(const string filename)const{
  OUTFILE(outf,filename.c_str());
  if(!outf.is_open()) return false;
  const bool succ = write(outf);
  outf.close();
  return succ;
}

bool PartialConstraints::write(ofstream &outf)const{
  
  outf << "FrameId " << _frameid << endl;
  outf << "ConstraintNodes "<< numConNodes() << endl;
  const vector<set<int> > &nodes = getConNodesSet();

  map<int,int> sortedNodes;
  int index = 0;
  BOOST_FOREACH(const set<int>& s, nodes){
	BOOST_FOREACH(const int nodeId, s){
	  sortedNodes.insert(make_pair<int,int>(nodeId,index));
	}
  }

  // BOOST_FOREACH(const map<int,int>::value_type &ele, sortedNodes){
  // 	const int nodeId = ele.first;
  // 	const int index = ele.second;
  // }

  // outf << i << " ";
  outf << "TargetPositions "<< _pc << endl;

  return outf.good();
}

bool PartialConstraints::load(const string filename){
  
  INFILE(inf,filename.c_str());
  if(!inf.is_open()) return false;
  const bool succ = load(inf);
  return succ;
}

bool PartialConstraints::load(ifstream &inf){
  
  string tempt_str;
  inf >> tempt_str >> _frameid;
  assert_ge(_frameid,0);
  
  int num_con_nodes = 0;
  inf >> tempt_str >> num_con_nodes;
  assert_ge(num_con_nodes,0);
  vector<int> nodes(num_con_nodes);
  for (int i = 0; i < num_con_nodes; ++i){
    inf >> nodes[i];
  }
  this->clear();
  addConNodes(nodes);
  
  inf >> tempt_str;
  _pc.resize(3,num_con_nodes);
  for (int r = 0; r < 3; ++r){
	for (int c = 0; c < num_con_nodes; ++c)
	  inf >> _pc(r,c);
  }
  
  return inf.good();
}

bool PartialConstraintsSet::load(const string filename){
  
  INFILE(inf,filename.c_str());
  if(!inf.is_open()) return false;

  string tempt_str;
  inf >> tempt_str;
  int number_of_frames=0;
  inf >> number_of_frames;
  bool succ = true;
  for (int i = 0; i < number_of_frames && succ; ++i){
    pPartialConstraints pcon = pPartialConstraints(new PartialConstraints());
	succ = pcon->load(inf);
	_partialCons.insert(pcon);
  }
  return succ;
}

bool PartialConstraintsSet::write(const string filename)const{
  
  OUTFILE(outf,filename.c_str());  
  if(!outf.is_open()) return false;

  outf << "frames " << _partialCons.size() << endl;
  bool succ = true;
  BOOST_FOREACH(pPartialConstraints pcon, _partialCons){
	succ = pcon->write(outf);
	if(!succ) break;
  }
  return succ;
}
