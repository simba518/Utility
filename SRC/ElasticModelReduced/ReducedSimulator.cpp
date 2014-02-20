#include "ReducedSimulator.h"
#include <JsonFilePaser.h>
#include <Log.h>
#include <ConMatrixTools.h>
#include <Timer.h>
using namespace UTILITY;
using namespace SIMULATOR;

bool ReducedSimulator::init(const string filename){

  JsonFilePaser jsonf;
  bool succ = false;
  if (jsonf.open(filename)){
	succ = true;
	jsonf.read("alpha_k",alpha_k,0.0);
	jsonf.read("alpha_m",alpha_m,0.0);
	jsonf.read("h",h,1.0);
  }else{
	ERROR_LOG("failed to open the initfile: " << filename);
  }
  return succ;
}

void ReducedImpLogConSimulator::setConGroups(const vector<int> &con_nodes){

  if (con_nodes.size() <= 0){
	removeAllCon();
	return ;
  }

  assert_ge(fullDim(), con_nodes.size()*3);
  SparseMatrix<double> sC;
  computeConM(con_nodes, sC, fullDim()/3);
  assert(model);
  assert_eq(sC.cols(), model->getModalBasis().rows());
  C = sC*model->getModalBasis();

  const int r = reducedDim();
  const int rc = C.rows();
  const MatrixXd hC = h*C;

  A.resize(r+rc,r+rc);
  A.setZero();
  A.topRightCorner(C.cols(),C.rows()) = hC.transpose();
  A.bottomLeftCorner(C.rows(),C.cols()) = hC;
  
  b.resize(r+rc);
}

void ReducedImpLogConSimulator::removeAllCon(){
  
  uc.resize(0);
  C.resize(0,0);
  const int r = reducedDim();
  A.resize(r,r);
  b.resize(r);
}

bool ReducedImpLogConSimulator::forward(){

  // FUNC_TIMER();

  assert(model);
  assert_eq(C.rows(), uc.size());

  static MatrixXd L;
  static VectorXd s;
  model->evaluateFK(q,s,L);

  L *= (h*h+h*alpha_k);
  L += (1.0f+h*alpha_m)*model->getReducedMassM();

  s *= (-1.0f*h);
  s += h*red_fext;
  s += model->getReducedMassM()*v;

  const int r = reducedDim();
  A.topLeftCorner(r,r) = L;
  b.head(r) = s;
  assert_eq(b.size(), r+C.rows());
  if (C.rows() > 0){
	b.tail(C.rows()) = uc-C*q;
  }

  v = A.lu().solve(b).head(r);
  q += h*v;

  return true;
}
