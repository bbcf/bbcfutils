// hit_annotator: scanning sequences with PWM
// Copyright (C) 2011, Eric R Paquet, Bernhard Sonderegger, Felix Naef, EPFL
#ifndef _FIXSIZESORTEDLIST_
#define _FIXSIZESORTEDLIST_

#include<list>
#include<limits>
#include<string>

using namespace std;

typedef pair<double,pair<int,string*> > seqScore;

const seqScore C_MAX(numeric_limits<double>::max(),make_pair(0,new string("")));
const seqScore C_MIN(numeric_limits<double>::min(),make_pair(0,new string("")));

class fsslist
{
 public:
  fsslist(int size):
    maxSize(size){
  }
  const seqScore& max()const{
    return elements.size() > 0 ? elements.back(): C_MAX;
  }
  const seqScore& min()const{
    return elements.size() > 0 ? elements.front(): C_MIN;
  }
  void insert(const seqScore& e){
    list<seqScore>::iterator it = elements.begin();
    for (;it != elements.end();++it){
      if (e.first < it->first)break;
    }
    itL.push_back(elements.insert(it,e));
    //it = elements.begin();
    //for (;it != elements.end();++it){
    //  cout << "["<<it->first << ","<<it->second.first<<"],";
    //}
    //cout<<endl;
    if (itL.size() > maxSize){
      delete itL.front()->second.second;
      elements.erase(itL.front());
      itL.pop_front();
    }
  }
  bool full()const{return elements.size()==maxSize;}
  void clear(){
    elements.clear();
    itL.clear();
  }
 private:
  list<seqScore> elements;
  list<list<seqScore>::iterator> itL;
  int maxSize;
};


#endif //_FIXSIZESORTEDLIST_
