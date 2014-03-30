#include <SparseGenEigenSolver.h>
#include <SparseMatrixTools.h>
#include <MassMatrix.h>
#include <ElasticForceTetFullStVK.h>
#include "nnls.h"
#include "MaterialFitting.h"
using namespace UTILITY;
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
  timer.stop("computeK: ");
}

void MaterialFitting::computeM(){
  
  TRACE_FUN();
  INFO_LOG("compute mass matrix");
  ComputeMassMat mass;
  initDensity(tetmesh->tets().size(),rho);
  Timer timer;
  timer.start();
  mass.compute(M,*tetmesh,rho);
  timer.stop("computeM: ");
}

void MaterialFitting::computeK(SparseMatrix<double> &K)const{
  
  ElasticForceTetFullStVK elas_c(tetmesh);
  const bool succ = elas_c.prepare(); assert(succ);
  VectorXd x0(tetmesh->nodes().size()*3);
  tetmesh->nodes(x0);
  K = elas_c.K(x0);
  K = -K;
}

void MaterialFitting::computeK(SparseMatrix<double> &K,const vector<double> &G,const vector<double> &L)const{

  const vector<double> G0 = tetmesh->material()._G;
  const vector<double> L0 = tetmesh->material()._lambda;
  
  tetmesh->material()._G = G;
  tetmesh->material()._lambda = L;

  ElasticForceTetFullStVK elas_c(tetmesh);
  const bool succ = elas_c.prepare(); assert(succ);
  VectorXd x0(tetmesh->nodes().size()*3);
  tetmesh->nodes(x0);
  K = elas_c.K(x0);
  K = -K;

  tetmesh->material()._G = G0;
  tetmesh->material()._lambda = L0;
}

void MaterialFitting::computeM(DiagonalMatrix<double,-1>&M)const{

  assert(tetmesh);
  MassMatrix mass;
  mass.compute(M,*tetmesh);
}

void MaterialFitting::computeM(DiagonalMatrix<double,-1>&M,const vector<double>&rh)const{
  
  assert(tetmesh);
  MassMatrix mass;
  const vector<double> r0 = tetmesh->material()._rho;
  tetmesh->material()._rho = rh;
  mass.compute(M,*tetmesh);
  tetmesh->material()._rho = r0;
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
  // addAverageDensityObjfun(objfun);
  // addFixedNodesObjfun(objfun);

  const SXMatrix M1 = assembleObjMatrix();
  for (int i = 0; i < M1.size1(); ++i){
    for (int j = 0; j < M1.size2(); ++j)
  	  objfun += M1.elem(i,j)*M1.elem(i,j);
  }
  objfun *= 0.5f;
}

void MaterialFitting::solveByNNLS(){

  const SXMatrix M1 = assembleObjMatrix();
  vector<SX> smooth_funs;
  computeSmoothObjFunctions(smooth_funs);
  
  SXMatrix obj_fun_m(M1.size1()*M1.size2()+smooth_funs.size(),1);
  for (int i = 0; i < M1.size1(); ++i){
    for (int j = 0; j < M1.size2(); ++j)
	  obj_fun_m.elem(i*M1.size2()+j,0) = M1.elem(i,j);
  }
  for (int i = 0; i < smooth_funs.size(); ++i){
    obj_fun_m.elem(M1.size2()*M1.size2()+i,0) = smooth_funs[i];
  }

  VSX variable_x;
  initAllVariables(variable_x);
  CasADi::SXFunction g_fun = CasADi::SXFunction(variable_x,obj_fun_m);
  g_fun.init();
  
  VectorXd x0(variable_x.size()), b;
  x0.setZero();
  CASADI::evaluate(g_fun, x0, b);
  b = -b;

  MatrixXd J = CASADI::convert<double>(g_fun.jac());
  assert_eq(J.rows(), b.size());
  assert_eq(J.cols(), x0.size());
  
  VectorXd x(b.size()),  w(J.cols()), zz(J.rows());
  VectorXi index(J.cols()*2);
  getInitValue(x);
  
  int exit_code = 0;
  double residual = 1e-18;
  nnls(&J(0,0), J.rows(), J.rows(), J.cols(), 
  	   &b[0], &x[0], &residual,
  	   &w[0], &zz[0], &index[0], &exit_code);
  INFO_LOG("residual: " << residual);
  print_NNLS_exit_code(exit_code);

  rlst.resize(x.size());
  for (int i = 0; i < x.size(); ++i){
    rlst[i] = x[i];
  }

  const MatrixXd H = J.transpose()*J;
  SelfAdjointEigenSolver<MatrixXd> es;
  es.compute(H);
  cout << "The eigenvalues of H are: " << es.eigenvalues().transpose() << endl;
}

void MaterialFitting::solveByIpopt(){
  
  TRACE_FUN();

  // init solver
  VSX variable_x;
  initAllVariables(variable_x);
  CasADi::SXFunction fun = CasADi::SXFunction(variable_x,objfun);
  CasADi::IpoptSolver solver = CasADi::IpoptSolver(fun);
  solver.setOption("generate_hessian",use_hessian);
  solver.setOption("tol",1e-18);
  solver.setOption("max_iter",20000);
  Timer timer;
  timer.start();
  solver.init();
  timer.stop("solver init: ");

  // set init value
  VectorXd init_x;
  getInitValue(init_x);
  solver.setInput(&init_x[0],CasADi::NLP_X_INIT);

  // set bounds
  vector<double> lower(variable_x.size()), upper(variable_x.size());
  for (int i = 0; i < lower.size(); ++i){
    lower[i] = lower_bound;
    upper[i] = upper_bound;
  }
  if (lower_bound >= 0.0f){
	solver.setInput(lower,CasADi::NLP_LBX);
  }
  if (upper_bound > 0.0f){
	assert_gt(upper_bound, lower_bound);
	solver.setInput(upper,CasADi::NLP_UBX);
  }

  // solving
  solver.solve();
  rlst.resize(variable_x.size());
  solver.getOutput(rlst,CasADi::NLP_X_OPT);
}

void MaterialFitting::solveByLinearSolver(){
  
  MatrixXd H;
  VectorXd g;
  hessGrad(H,g);

  SelfAdjointEigenSolver<MatrixXd> es;
  es.compute(H);
  cout << "H.norm() = " << H.norm() << endl;
  cout << "The eigenvalues of H are: \n" << es.eigenvalues().transpose() << endl;

  const VectorXd x = H.lu().solve(-g);
  rlst.resize(x.size());
  for (int i = 0; i < x.size(); ++i)
    rlst[i] = x[i];
  cout<< "\n(H*x+g).norm() = " << (H*x+g).norm() << "\n\n";
}

SXMatrix MaterialFitting::assembleObjMatrix(){
  const SXMatrix M1 = (K.mul(W)-M.mul(W.mul(lambda)));
  return M1;
}

void MaterialFitting::computeSmoothObjFunctions(vector<SX> &funs)const{

  const VVec4i &neigh_tet = tetmesh->faceNeighTet();
  funs.clear();
  funs.reserve(neigh_tet.size()*4);

  const VSX G = getG();
  const VSX Lame = getLame();
  const VSX rho = getRho();
  const VSX E = getE();
  const VSX v = getV();
  for (int i = 0; i < neigh_tet.size(); ++i){

	const double vi = tetmesh->volume(i);
	assert_gt(vi,0);
	for (int k = 0; k < 4; ++k){
	  const int j = neigh_tet[i][k];
	  if (j >= 0){
		const double vj = tetmesh->volume(j);
		assert_gt(vj,0);
		assert(i!=j);
		if (mu_neigh_G > 0)
		  funs.push_back( 0.5f*sqrt(mu_neigh_G*(vi+vj))*(G[i]-G[j]) );
		if (mu_neigh_L > 0)
		  funs.push_back( 0.5f*sqrt(mu_neigh_L*(vi+vj))*(Lame[i]-Lame[j]) );
		if (mu_neigh_rho > 0)
		  funs.push_back( 0.5f*sqrt(mu_neigh_rho*(vi+vj))*(rho[i]-rho[j]) );
		if (mu_neigh_E > 0)
		  funs.push_back( 0.5f*sqrt(mu_neigh_E*(vi+vj))*(E[i]-E[j]) );
		if (mu_neigh_poission > 0)
		  funs.push_back( 0.5f*sqrt(mu_neigh_poission*(vi+vj))*(v[i]-v[j]) );
	  }
	}
  }
}

void MaterialFitting::addSmoothObjfun(SX &objfun)const{
  
  vector<SX> funs;
  computeSmoothObjFunctions(funs);
  for (int i = 0; i < funs.size(); ++i){
    objfun += funs[i]*funs[i];
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
		  const double scalor = 1e-8; /// @todo how to choose scalor?
		  objfun += 0.25f*(vi+vj)*(G[i]-G[j])*(G[i]-G[j])*scalor;
		  objfun += 0.25f*(vi+vj)*(Lame[i]-Lame[j])*(Lame[i]-Lame[j])*scalor;
		  objfun += 0.25f*(vi+vj)*(rho[i]-rho[j])*(rho[i]-rho[j])*scalor;
		}
	  }
	}
  }
}

void MaterialFitting::hessGrad(MatrixXd &H, VectorXd &g)const{
  
  VSX variable_x;
  initAllVariables(variable_x);
  CasADi::SXFunction fun = CasADi::SXFunction(variable_x,objfun);
  fun.init();

  CasADi::SXFunction grad_fun = CasADi::SXFunction(variable_x,fun.jac());
  grad_fun.init();

  CasADi::SXFunction hess_fun = CasADi::SXFunction(variable_x,fun.hess());
  hess_fun.init();

  VectorXd x0(variable_x.size());
  x0.setZero();
  CASADI::evaluate(grad_fun, x0, g);
  CASADI::evaluate(hess_fun, x0, H);

  assert_eq(g.norm(), g.norm());
  assert_eq(H.norm(), H.norm());
}

void MaterialFitting::print_NNLS_exit_code(const int exit_code)const{
  
  switch (exit_code){
  case 1: INFO_LOG("nnls solve success!");
	break;
  case 2: ERROR_LOG("nnls failed, the dimensions of the problem are bad, either m .le. 0 or n .le. 0.");
	break;
  case 3: ERROR_LOG("nnls failed, iteration count exceeded.");
	break;
  default: ERROR_LOG("nnls failed, unknown reason.");
	break;
  }
}

void MaterialFitting::testFixedNodes(){
  
  if(fixednodes.size() <= 0){
	return ;
  }

  const VVec4i &neigh_tet = tetmesh->faceNeighTet();
  const VVec4i &tets = tetmesh->tets();
  for (int i = 0; i < tets.size(); ++i){
	int findnodes = 0;
    for (int j = 0; j < 4; ++j){
	  if ( fixednodes.find(tets[i][j])!=fixednodes.end() )
		findnodes ++;
	}
	if (4 == findnodes){
	  tetmesh->material()._rho[i] *= 100.0f;
	  tetmesh->material()._G[i] *= 100.0f;
	  tetmesh->material()._lambda[i] *= 100.0f;
	}
  }
}

void MaterialFitting::saveResults(const string filename)const{
  
  const ElasticMaterial<double> mtl_0 = tetmesh->material();
  const vector<double> g = getShearGResult();
  const vector<double> l = getLameResult();
  const vector<double> r = getDensityResult();

  tetmesh->material()._G = g;
  tetmesh->material()._lambda = l;
  tetmesh->material()._rho = r;

  const bool succ = tetmesh->writeElasticMtlVTK(filename); assert(succ);

  tetmesh->material() = mtl_0;
}

void MaterialFitting::printResult()const{
  
  cout << "\nnumber of variables: " << rlst.size() << endl;
  for (int i = 0; i < rlst.size(); ++i){
    cout << rlst[i]<< " ";
  }
  cout << "\n\n";
}

void MaterialFitting_MA_K::initDensity(const int num_tet, VSX &rho)const{

  const vector<double> &r = tetmesh->material()._rho;
  assert_eq(num_tet, r.size());
  rho.resize(num_tet);
  for (int i = 0; i < num_tet; ++i)	rho[i] = r[i];
}

void MaterialFitting_MA_K::getInitValue(VectorXd &init_x)const{

  const vector<double> &g = tetmesh->material()._G;
  const vector<double> &la = tetmesh->material()._lambda;
  const int num_tet = tetmesh->tets().size();
  init_x.resize(num_tet*2);
  for (int i = 0; i < num_tet; ++i){
    init_x[i] = g[i];
    init_x[i+num_tet] = la[i];
  }
}

void MaterialFitting_Diag_KM::assembleObjfun(){
  
  TRACE_FUN();
  objfun = 0.0f;
  addSmoothObjfun(objfun);
  addAverageDensityObjfun(objfun);
  addFixedNodesObjfun(objfun);

  const MatrixXd I = MatrixXd::Identity(lambda.size(),lambda.size());
  const SXMatrix M1 = trans(W).mul(K.mul(W))-lambda;
  const SXMatrix _Diag_KM = trans(W).mul(M.mul(W))-CASADI::convert(I);
  for (int i = 0; i < M1.size1(); ++i){
    for (int j = 0; j < M1.size2(); ++j){
	  objfun += mu_stiff*M1.elem(i,j)*M1.elem(i,j);
	  objfun += mu_mass*_Diag_KM.elem(i,j)*_Diag_KM.elem(i,j);
	}
  }

  objfun *= 0.5f;
}

SXMatrix MaterialFitting_Diag_M::assembleObjMatrix(){

  const MatrixXd I = MatrixXd::Identity(lambda.size(),lambda.size());
  const SXMatrix M2 = trans(W).mul(M.mul(W))-CASADI::convert(I);
  return M2;
}

void MaterialFitting_Diag_M::initShearG(const int num_tet, VSX &G){
  
  const vector<double> &r = tetmesh->material()._G;
  assert_eq(num_tet, r.size());
  G.resize(num_tet);
  for (int i = 0; i < num_tet; ++i){
	G[i] = r[i];
  }
}

void MaterialFitting_Diag_M::initLame(const int num_tet, VSX &Lame){

  const vector<double> &r = tetmesh->material()._lambda;
  assert_eq(num_tet, r.size());
  Lame.resize(num_tet);
  for (int i = 0; i < num_tet; ++i){
	Lame[i] = r[i];
  }
}

void MaterialFitting_Diag_M::getInitValue(VectorXd &init_x)const{
    
  const vector<double> &r = tetmesh->material()._rho;
  const int num_tet = tetmesh->tets().size();
  init_x.resize(num_tet);
  for (int i = 0; i < num_tet; ++i){
    init_x[i] = r[i];
  }
}

vector<double> MaterialFitting_Diag_M::getDensityResult()const{

  const int num_tet = tetmesh->tets().size();
  assert_ge(rlst.size(),num_tet);
  return vector<double>(rlst.begin(),rlst.begin()+num_tet);
}

SXMatrix MaterialFitting_Diag_K::assembleObjMatrix(){

  const SXMatrix M1 = trans(W).mul(K.mul(W))-lambda;
  return M1;
}

void MaterialFitting_EV::getInitValue(VectorXd &init_x)const{
  
  const vector<double> &g = tetmesh->material()._G;
  const vector<double> &la = tetmesh->material()._lambda;
  const vector<double> E = MaterialFitting::getYoungE(g,la);
  const vector<double> v = MaterialFitting::getPoissonV(g,la);
  
  const int num_tet = tetmesh->tets().size();
  init_x.resize(num_tet*2);
  for (int i = 0; i < num_tet; ++i){
    init_x[i] = E[i];
    init_x[i+num_tet] = v[i];
  }
}

SXMatrix MaterialFitting_EV::assembleObjMatrix(){

  const SXMatrix M1 = trans(W).mul(K.mul(W))-lambda;
  return M1;
}
