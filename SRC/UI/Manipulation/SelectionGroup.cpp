#include <boost/foreach.hpp>
#include "SelectionGroup.h"
using namespace UTILITY;

void SelectionGroup::addGroup( const set<int>& g){

  if(!g.empty()){
	removeGroup(g);
	groups.push_back(g);
  }
}

void SelectionGroup::addGroup( const vector<int>& g){

  addGroup(vec2set(g));
}

void SelectionGroup::removeGroup( const set<int>& remove_node_ids){

  vector<set<int> > tempt_con_group;

  vector<set<int> >::iterator g = groups.begin();
  for(; g!=groups.end(); g++ ){

	set<int>::const_iterator r = remove_node_ids.begin();
	for ( ;r!= remove_node_ids.end(); r++){
	  g->erase(*r);
	}
	if(!g->empty()){
	  tempt_con_group.push_back(*g);
	}
  }
  this->groups = tempt_con_group; 
}

void SelectionGroup::removeGroup( const vector<int>& g){

  removeGroup(vec2set(g));
}

void SelectionGroup::setGroup(const vector<set<int> > &groups){

  BOOST_FOREACH(const set<int> &ele, groups){
	addGroup(ele);
  }
}

void SelectionGroup::setGroup(const vector<vector<int> > &groups){

  BOOST_FOREACH(const vector<int> &ele, groups){
	addGroup(ele);
  }  
}

vector<vector<int> > SelectionGroup::getGroupVec()const{
  
  vector<vector<int> > g;
  BOOST_FOREACH(set<int> ele, groups){
	g.push_back(set2vec(ele));
  }
  return g;
}

set<int> SelectionGroup::vec2set(const vector<int> &v)const{
  
  set<int> s;
  BOOST_FOREACH(int ele, v){
	s.insert(ele);
  }
  return s;
}

vector<int> SelectionGroup::set2vec(const set<int> &s)const{

  vector<int> v;
  BOOST_FOREACH(int ele, s){
	v.push_back(ele);
  }
  return v;
}

bool SelectionGroup::contain(const int ele_id)const{
  
  BOOST_FOREACH(const set<int>& g, groups){
	if (g.find(ele_id) != g.end()){
	  return true;
	}
  }  
  return false;
}
