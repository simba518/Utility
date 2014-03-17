#include <SparseGenEigenSolver.h>
#include <SparseMatrixTools.h>
#include "MaterialFitting.h"
using namespace ELASTIC_OPT;

void MaterialFitting::computeK(){
  
  TRACE_FUN();
  ComputeStiffnessMat elas(tetmesh);
  initShearG(tetmesh->tets().size(),G);
  initLame(tetmesh->tets().size(),Lame);
  elas.setMaterial(G, Lame);
  bool succ = elas.prepare();  assert(succ);
  VectorXd x0;
  tetmesh->nodes(x0);

  Timer timer;
  timer.start();
  K = elas.K(x0);
  K *= -1.0f;
  timer.stop("elas.K(x0) ");
}

void MaterialFitting::computeM(){
  
  TRACE_FUN();
  INFO_LOG("compute mass matrix");
  ComputeMassMat mass;
  initDensity(tetmesh->tets().size(),rho);
  mass.compute(M,*tetmesh,rho);
}

void MaterialFitting::removeFixedDOFs(){
  
  if (fixednodes.size() > 0){

	SparseMatrix<double> Pm;
	EIGEN3EXT::genReshapeMatrix(K.size1(),3,fixednodes,Pm);
	SXMatrix P;
	CASADI::convert(Pm,P);
	assert_eq(P.size1(),K.size1()-fixednodes.size()*3);
	assert_eq(P.size2(),K.size1());
	K = P.mul(K.mul(trans(P)));
	M = P.mul(M.mul(trans(P)));
	W = P.mul(W);
  }
}

void MaterialFitting::getInitValue(VectorXd &init_x)const{
  
  const vector<double> &g = tetmesh->material()._G;
  const vector<double> &la = tetmesh->material()._lambda;
  const vector<double> &rh = tetmesh->material()._rho;
  const int num_tet = tetmesh->tets().size();

  init_x.resize(num_tet*3);
  for (int i = 0; i < num_tet; ++i){
    init_x[i] = g[i];
    init_x[i+num_tet] = la[i];
    init_x[i+num_tet*2] = rh[i];
  }
}

vector<double> MaterialFitting::getShearGResult()const{

  const int num_tet = tetmesh->tets().size();
  assert_ge(rlst.size(),num_tet);
  return vector<double>(rlst.begin(),rlst.begin()+num_tet);
}

vector<double> MaterialFitting::getLameResult()const{

  const int num_tet = tetmesh->tets().size();
  assert_ge(rlst.size(),num_tet*2);
  return vector<double>(rlst.begin()+num_tet,rlst.begin()+2*num_tet);
}

vector<double> MaterialFitting::getDensityResult()const{

  const int num_tet = tetmesh->tets().size();
  assert_ge(rlst.size(),num_tet*3);
  return vector<double>(rlst.begin()+2*num_tet,rlst.begin()+3*num_tet);
}

void MaterialFitting::assembleObjfun(){
  
  TRACE_FUN();
  objfun = 0.0f;
  addSmoothObjfun(objfun);
  addAverageDensityObjfun(objfun);
  addFixedNodesObjfun(objfun);

  const SXMatrix M1 = K.mul(W)-M.mul(W.mul(lambda));
  for (int i = 0; i < M1.size1(); ++i){
    for (int j = 0; j < M1.size2(); ++j)
	  objfun += M1.elem(i,j)*M1.elem(i,j);
  }

  objfun *= 0.5f;
}

void MaterialFitting::addSmoothObjfun(SX &objfun)const{
  
  const VVec4i &neigh_tet = tetmesh->faceNeighTet();
  for (int i = 0; i < neigh_tet.size(); ++i){
	const double vi = tetmesh->volume(i);
	assert_gt(vi,0);
	for (int k = 0; k < 4; ++k){
	  const int j = neigh_tet[i][k];
	  if (j >= 0){
		const double vj = tetmesh->volume(j);
		assert_gt(vj,0);
		assert(i!=j);
		objfun += 0.25f*mu_neigh_G*(vi+vj)*(G[i]-G[j])*(G[i]-G[j]);
		objfun += 0.25f*mu_neigh_L*(vi+vj)*(Lame[i]-Lame[j])*(Lame[i]-Lame[j]);
		objfun += 0.25f*mu_neigh_rho*(vi+vj)*(rho[i]-rho[j])*(rho[i]-rho[j]);
	  }
	}
  }
}

void MaterialFitting::addAverageDensityObjfun(SX &objfun)const{
  
  const vector<double> &desired_rho = tetmesh->material()._rho;
  const int num_tet = tetmesh->tets().size();
  assert_eq(desired_rho.size(), num_tet);
  assert_eq(rho.size(), num_tet);
  double desired_mass = 0.0f;
  SX mass = 0.0f;
  for (int i = 0; i < num_tet; ++i){
    const double v = tetmesh->volume(i);
	assert_gt(v,0.0f);
	desired_mass += desired_rho[i]*v;
	mass += rho[i]*v;
  }
  assert_gt(desired_mass, 0.0f);
  
  objfun += mu_average_rho*(desired_mass-mass)*(desired_mass-mass);
}

void MaterialFitting::addFixedNodesObjfun(SX &objfun)const{
  
  if(fixednodes.size() <= 0){
	return ;
  }

  const VVec4i &neigh_tet = tetmesh->faceNeighTet();
  const VVec4i &tets = tetmesh->tets();
  assert_eq(neigh_tet.size(), tets.size());
  assert_eq(G.size(), tets.size());
  assert_eq(Lame.size(), tets.size());
  assert_eq(rho.size(), tets.size());
  for (int i = 0; i < tets.size(); ++i){
	int findnodes = 0;
    for (int j = 0; j < 4; ++j){
	  if ( fixednodes.find(tets[i][j])!=fixednodes.end() ){
		findnodes ++;
	  }
	}
	if (4 == findnodes){
	  const double vi = tetmesh->volume(i);
	  for (int k = 0; k < 4; ++k){
		const int j = neigh_tet[i][k];
		if (j >= 0){
		  const double vj = tetmesh->volume(j);
		  assert_gt(vj,0);
		  assert(i!=j);
		  const double scalor = 1e-16; /// @todo how to choose scalor?
		  objfun += 0.25f*(vi+vj)*(G[i]-G[j])*(G[i]-G[j])*scalor;
		  objfun += 0.25f*(vi+vj)*(Lame[i]-Lame[j])*(Lame[i]-Lame[j])*scalor;
		  objfun += 0.25f*(vi+vj)*(rho[i]-rho[j])*(rho[i]-rho[j])*scalor;
		}
	  }
	}
  }
}

void MaterialFitting::solve(){
  
  TRACE_FUN();

  // init solver
  VSX variable_x;
  initAllVariables(variable_x);
  CasADi::SXFunction fun = CasADi::SXFunction(variable_x,objfun);
  CasADi::IpoptSolver solver = CasADi::IpoptSolver(fun);
  solver.setOption("generate_hessian",use_hessian);
  solver.setOption("tol",1e-18);
  solver.setOption("max_iter",1000);
  solver.init();

  // set init value
  VectorXd init_x;
  getInitValue(init_x);
  solver.setInput(&init_x[0],CasADi::NLP_X_INIT);

  // set bounds
  vector<double> lower;
  lower.resize(variable_x.size());
  for (int i = 0; i < lower.size(); ++i)
    lower[i] = 0.0f;
  solver.setInput(lower,CasADi::NLP_LBX);

  // solving
  solver.solve();
  rlst.resize(variable_x.size());
  solver.getOutput(rlst,CasADi::NLP_X_OPT);
}

void MaterialFitting::saveResults(const string filename)const{
  
  const ElasticMaterial<double> mtl_0 = tetmesh->material();

  tetmesh->material()._G = getShearGResult();
  tetmesh->material()._lambda = getLameResult();
  tetmesh->material()._rho = getDensityResult();

  const bool succ = tetmesh->writeElasticMtlVTK(filename); assert(succ);
  tetmesh->material() = mtl_0;
}

void MaterialFittingM2::assembleObjfun(){
  
  TRACE_FUN();
  objfun = 0.0f;
  addSmoothObjfun(objfun);
  addAverageDensityObjfun(objfun);
  addFixedNodesObjfun(objfun);

  const MatrixXd I = MatrixXd::Identity(lambda.size(),lambda.size());
  const SXMatrix M1 = trans(W).mul(K.mul(W))-lambda;
  const SXMatrix M2 = trans(W).mul(M.mul(W))-CASADI::convert(I);
  for (int i = 0; i < M1.size1(); ++i){
    for (int j = 0; j < M1.size2(); ++j){
	  objfun += mu_stiff*M1.elem(i,j)*M1.elem(i,j);
	  objfun += mu_mass*M2.elem(i,j)*M2.elem(i,j);
	}
  }

  objfun *= 0.5f;
}

void MaterialFittingM3::initDensity(const int num_tet, VSX &rho)const{

  const vector<double> &r = tetmesh->material()._rho;
  assert_eq(num_tet, r.size());
  rho.resize(num_tet);
  for (int i = 0; i < num_tet; ++i){
	rho[i] = r[i];
  }
}

void MaterialFittingM3::getInitValue(VectorXd &init_x)const{

  const vector<double> &g = tetmesh->material()._G;
  const vector<double> &la = tetmesh->material()._lambda;
  const int num_tet = tetmesh->tets().size();
  init_x.resize(num_tet*2);
  for (int i = 0; i < num_tet; ++i){
    init_x[i] = g[i];
    init_x[i+num_tet] = la[i];
  }
}
