#include <iostream>
#include <MatrixTools.h>
#include <MatrixIO.h>
#include <interfaces/ipopt/ipopt_solver.hpp>
#include <CASADITools.h>
#include <MassMatrix.h>
#include <SparseGenEigenSolver.h>
#include <SparseMatrixTools.h>
#include <ElasticForceTetFullStVK.h>
#include <Timer.h>
#include <MassMatrix.h>
#include "ComputeStiffnessMat.h"
#include "ComputeMassMat.h"
using namespace std;
using namespace UTILITY;
using namespace EIGEN3EXT;
using namespace ELASTIC_OPT;
using CasADi::SX;

void testK(){
  
  const string tet_fname = std::string(TEST_DATA_DIR)+"beam.abq";
  pTetMesh mesh = pTetMesh(new TetMesh);
  assert ( mesh->load(tet_fname) );
  mesh->material().reset(1000.0f,2E6,0.45);
  const int n = mesh->nodes().size();

  const string stiff_mat_f = std::string(TEST_DATA_DIR)+"beam_sparse_K.b";
  SparseMatrix<double> ref_K;
  assert ( load(ref_K,stiff_mat_f) );
  assert_eq(ref_K.rows(),n*3);
  assert_eq(ref_K.cols(),n*3);

  ComputeStiffnessMat elas(mesh);
  assert(elas.prepare());
  VectorXd x0(n*3);
  mesh->nodes(x0);
  const SXMatrix &K = elas.K(x0);
  
  MatrixXd Km;
  CASADI::convert(K,Km);
  const MatrixXd ref_Km = ref_K;
  
  assert_eq(Km.rows(),n*3);
  assert_eq(Km.cols(),n*3);
  assert_lt((Km+ref_Km).norm(),1e-5);
  
  cout<< "(K-K_ref).norm(): " << (Km+ref_Km).norm() << endl;
  cout << "test end.." << endl;
}

void testM(){
  
  const string tet_fname = std::string(TEST_DATA_DIR)+"two_tet.abq";
  const string tet_mtl_fname = std::string(TEST_DATA_DIR)+"two_tet.elastic";
  pTetMesh mesh = pTetMesh(new TetMesh);
  assert ( mesh->load(tet_fname) );
  assert ( mesh->loadElasticMtl(tet_mtl_fname) );
  const int n = mesh->nodes().size();

  MassMatrix mass_c;
  Eigen::DiagonalMatrix<double,-1> diagM;
  mass_c.compute(diagM, *mesh);

  ComputeMassMat mass;
  SXMatrix M;
  mass.compute(M, *mesh);
  MatrixXd Me;
  CASADI::convert(M,Me);

  assert_eq(Me.rows(), Me.cols());
  assert_eq(Me.rows(), diagM.rows());
  assert_eq(Me.rows(), diagM.cols());
  const MatrixXd mm = Me-MatrixXd(diagM);
  assert_lt(mm.norm(),1e-8);
  cout<< "(M-Mref).norm() = " << mm.norm() << endl;
}

void testKSMall(){
  
  const string tet_fname = std::string(TEST_DATA_DIR)+"one_tet.abq";
  const string tet_mtl_fname = std::string(TEST_DATA_DIR)+"one_tet.elastic";
  pTetMesh mesh = pTetMesh(new TetMesh);
  assert ( mesh->load(tet_fname) );
  assert ( mesh->loadElasticMtl(tet_mtl_fname) );
  const int n = mesh->nodes().size();

  VectorXd x0(n*3);
  mesh->nodes(x0);

  ElasticForceTetFullStVK elas_c(mesh);
  assert(elas_c.prepare());
  const MatrixXd ref_Km = elas_c.K(x0);

  ComputeStiffnessMat elas(mesh);
  assert(elas.prepare());
  elas.K(x0);
  const SXMatrix &K = elas.K(x0);
  MatrixXd Km;
  CASADI::convert(K,Km);
  
  assert_eq(Km.rows(),n*3);
  assert_eq(Km.cols(),n*3);
  assert_lt((Km-ref_Km).norm(),1e-5);
  
  cout<< "(K-K_ref).norm(): " << (Km-ref_Km).norm() << endl;
  cout << "test end.." << endl;
}

void computeEigenValues(const string data_root,const int eigenNum,const set<int> &fixednodes){
  
  // load data.
  INFO_LOG("load data");
  pTetMesh tetmesh = pTetMesh(new TetMesh());
  bool succ = tetmesh->load(data_root+"mesh.abq"); assert(succ);
  succ = tetmesh->loadElasticMtl(data_root+"mesh.elastic"); assert(succ);
  const int num_tet = (int)tetmesh->tets().size();
  const int n = tetmesh->nodes().size();

  // save material
  INFO_LOG("save material");
  succ = tetmesh->writeElasticMtlVTK("./tempt/material_correct"); assert(succ);

  // compute K
  INFO_LOG("compute K");
  ElasticForceTetFullStVK elas(tetmesh);
  assert(elas.prepare());
  VectorXd x0(n*3);
  tetmesh->nodes(x0);
  SparseMatrix<double> K = elas.K(x0);
  K *= -1.0f;

  // compute mass matrix
  INFO_LOG("compute mass matrix");
  MassMatrix mass;
  DiagonalMatrix<double,-1> M;
  mass.compute(M,*tetmesh);

  // remove fixed nodes
  SparseMatrix<double> P;
  EIGEN3EXT::genReshapeMatrix(K.rows(),3,fixednodes,P);
  K = P*(K*P.transpose());
  const SparseMatrix<double> diagM = P*(M*P.transpose());

  /// compute W, lambda
  INFO_LOG("compute W, lambda");
  MatrixXd W_full;
  VectorXd lambda_full;
  const SparseMatrix<double> Klower = EIGEN3EXT::getLower(K);
  succ = EigenSparseGenEigenSolver::solve(Klower,diagM,W_full,lambda_full,eigenNum); 
  assert(succ);

  // const MatrixXd W = W_full.rightCols(lambda_full.size()-6);
  // const VectorXd lambda = lambda_full.tail(lambda_full.size()-6);

  // test Wt*M*W
  INFO_LOG("check results");
  const MatrixXd WtMWI = W_full.transpose()*diagM*W_full-MatrixXd::Identity(W_full.cols(),W_full.cols());
  const MatrixXd KW_MWLambda = K*W_full-(diagM*W_full)*lambda_full.asDiagonal();
  cout << "(WtMW-I).norm(): " << WtMWI.norm() << endl;
  cout << "(KW_MWLambda).norm(): " << KW_MWLambda.norm() << endl;
  cout << "eigenvalues: " << lambda_full.transpose() << endl;
  cout << "norm(Klower): " << Klower.norm() << endl;
  cout << "norm(M): " << diagM.diagonal().norm() << endl;

  // add fixed nodes to W, and lambda, then save
  const MatrixXd W = P.transpose()*W_full;
  const VectorXd lambda = lambda_full;
  succ = write(data_root+"tempt_eigenvalues.b", lambda); assert(succ);
  succ = write(data_root+"tempt_eigenvectors.b", W); assert(succ);
  
}

void recoverFullMtl(pTetMesh tetmesh, const set<int> &fixednodes,
					const MatrixXd &eigen_W, const VectorXd &eigen_lambda,
					const double mu_neigh, const string save_to){

  // load data.
  INFO_LOG("loading data");
  const int num_tet = (int)tetmesh->tets().size();
  const int n = tetmesh->nodes().size();

  SXMatrix W, lambda;
  CASADI::convert(eigen_W, W);
  assert_eq(W.size1(), n*3);
  assert_lt(W.size2(), W.size1());

  lambda.resize(eigen_lambda.size(), eigen_lambda.size());
  for (int i = 0; i < eigen_lambda.size(); ++i){
    lambda.elem(i,i) = eigen_lambda[i];
  }
  assert_eq(lambda.size1(), W.size2());
  assert_eq(lambda.size2(), W.size2());
  
  // init variables
  INFO_LOG("init variables");
  const CASADI::VSX G = CASADI::makeSymbolic(num_tet, "G");
  const CASADI::VSX Lame = CASADI::makeSymbolic(num_tet, "lame");
  const CASADI::VSX G_Lame = CASADI::connect(G,Lame);
  assert_eq(G_Lame.size(), num_tet*2);

  // compute K
  INFO_LOG("compute K");
  ComputeStiffnessMat elas(tetmesh);
  elas.setMaterial(G, Lame);
  assert(elas.prepare());
  VectorXd x0(n*3);
  tetmesh->nodes(x0);

  Timer timer;
  timer.start();
  SXMatrix K = elas.K(x0);
  timer.stop("elas.K(x0) ");

  assert_eq(K.size1(), n*3);
  assert_eq(K.size2(), n*3);
  K *= -1.0f;

  // cout << "K:\n";
  // cout << K << endl;
  // cout << endl << endl;

  // compute mass matrix
  INFO_LOG("compute mass matrix");
  MassMatrix mass;
  DiagonalMatrix<double,-1> diag_M;
  mass.compute(diag_M,*tetmesh);
  SXMatrix M(diag_M.rows(), diag_M.cols());
  for (int i = 0; i < diag_M.rows(); ++i){
    M.elem(i,i) = diag_M.diagonal()[i];
  }
  assert_eq(M.size1(), n*3);
  assert_eq(M.size2(), n*3);


  /// remove fixed nodes
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

  // assemble objective function.
  INFO_LOG("assemble objective function.");
  SX objfun = 0.0f;
  const SXMatrix M1 = K.mul(W)-M.mul(W.mul(lambda));
  for (int i = 0; i < M1.size1(); ++i){
    for (int j = 0; j < M1.size2(); ++j)
	 objfun += M1.elem(i,j)*M1.elem(i,j);
  }

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
		objfun += 0.5f*mu_neigh*(vi+vj)*(G[i]-G[j])*(G[i]-G[j]);
		objfun += 0.5f*mu_neigh*(vi+vj)*(Lame[i]-Lame[j])*(Lame[i]-Lame[j]);
	  }
	}
  }
  objfun *= 0.5f;

  // solve.
  INFO_LOG("init solver");
  CasADi::SXFunction fun = CasADi::SXFunction(G_Lame,objfun);
  CasADi::IpoptSolver solver = CasADi::IpoptSolver(fun);

  solver.setOption("generate_hessian",true);
  solver.setOption("tol",1e-18);
  solver.setOption("max_iter",1000);

  timer.start();
  solver.init();
  timer.stop("solver.init() ");

  VectorXd init_x;
  init_x.resize(G_Lame.size());
  const vector<double> &g = tetmesh->material()._G;
  const vector<double> &la = tetmesh->material()._lambda;
  for (int i = 0; i < num_tet; ++i){
    init_x[i] = g[i];
    init_x[i+num_tet] = la[i];
  }
  solver.setInput(&init_x[0],CasADi::NLP_X_INIT);

  cout << endl << endl;
  INFO_LOG("initial x: ");
  for (int i = 0; i < init_x.size(); ++i){
    cout << init_x[i] << " ,";
  }
  cout << endl << endl;

  vector<double> lower;
  lower.resize(G_Lame.size());
  for (int i = 0; i < num_tet; ++i){
    lower[i] = 0.0f;
    lower[i+num_tet] = 0.0f;
  }
  solver.setInput(lower,CasADi::NLP_LBX);

  // solving
  INFO_LOG("solving");
  timer.start();
  solver.solve();
  timer.stop("solver.solve() ");
  std::vector<double> vx(G_Lame.size());
  solver.getOutput(vx,CasADi::NLP_X_OPT);  

  // convert to Young's(E) and Poisson(v)
  vector<double> Young_E(num_tet), Poisson_v(num_tet);
  for (int i = 0; i < num_tet; ++i){

	const double gi = vx[i];
	const double li = vx[i+num_tet];
	const Matrix<double,2,1> Ev = ElasticMaterial<double>::fromLameConstant(gi,li);
	Young_E[i] = Ev(0,0);
	Poisson_v[i] = Ev(1,0);
	
	tetmesh->material()._G[i] = gi;
	tetmesh->material()._lambda[i] = li;
  }

  // save results
  INFO_LOG("save results");

  cout << endl << endl << "G, L: \n";
  for (int i = 0; i < num_tet; ++i){
	cout << vx[i] << ", " << vx[i+num_tet] << endl;
  }
  cout << endl << endl;

  cout << endl << endl << "E, v: \n";
  for (int i = 0; i < num_tet; ++i){
	cout << Young_E[i] << ", " << Poisson_v[i] << endl;
  }
  cout << endl << endl;

  bool succ = tetmesh->writeElasticMtlVTK(save_to);
  assert(succ);
}

void recoverFullMtlMethod2(pTetMesh tetmesh, const set<int> &fixednodes,
						   const MatrixXd &eigen_W, const VectorXd &eigen_lambda,
						   const double mu_neigh, const string save_to){

  // load data.
  INFO_LOG("loading data");
  const int num_tet = (int)tetmesh->tets().size();
  const int n = tetmesh->nodes().size();

  SXMatrix W, lambda;
  CASADI::convert(eigen_W, W);
  assert_eq(W.size1(), n*3);
  assert_lt(W.size2(), W.size1());

  lambda.resize(eigen_lambda.size(), eigen_lambda.size());
  for (int i = 0; i < eigen_lambda.size(); ++i){
    lambda.elem(i,i) = eigen_lambda[i];
  }
  assert_eq(lambda.size1(), W.size2());
  assert_eq(lambda.size2(), W.size2());
  
  // init variables
  INFO_LOG("init variables");
  const CASADI::VSX G = CASADI::makeSymbolic(num_tet, "G");
  const CASADI::VSX Lame = CASADI::makeSymbolic(num_tet, "lame");
  const CASADI::VSX rho = CASADI::makeSymbolic(num_tet, "rho");
  const CASADI::VSX G_Lame_rho = CASADI::connect(CASADI::connect(G,Lame),rho);

  // compute K
  INFO_LOG("compute K");
  ComputeStiffnessMat elas(tetmesh);
  elas.setMaterial(G, Lame);
  assert(elas.prepare());
  VectorXd x0(n*3);
  tetmesh->nodes(x0);

  Timer timer;
  timer.start();
  SXMatrix K = elas.K(x0);
  K *= -1.0f;
  timer.stop("elas.K(x0) ");

  // compute mass matrix
  INFO_LOG("compute mass matrix");
  ComputeMassMat mass;
  SXMatrix M;
  mass.compute(M,*tetmesh,rho);

  /// remove fixed nodes
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

  // assemble objective function.
  INFO_LOG("assemble objective function.");

  SX objfun = 0.0f;
  const SXMatrix M1 = K.mul(W)-M.mul(W.mul(lambda));
  for (int i = 0; i < M1.size1(); ++i){
    for (int j = 0; j < M1.size2(); ++j)
	  objfun += M1.elem(i,j)*M1.elem(i,j);
  }

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
		objfun += 0.25f*mu_neigh*(vi+vj)*(G[i]-G[j])*(G[i]-G[j]);
		objfun += 0.25f*mu_neigh*(vi+vj)*(Lame[i]-Lame[j])*(Lame[i]-Lame[j]);
		objfun += 0.25f*mu_neigh*(vi+vj)*(rho[i]-rho[j])*(rho[i]-rho[j]);
	  }
	}
  }
  objfun *= 0.5f;

  // solve.
  INFO_LOG("init solver");
  CasADi::SXFunction fun = CasADi::SXFunction(G_Lame_rho,objfun);
  CasADi::IpoptSolver solver = CasADi::IpoptSolver(fun);

  solver.setOption("generate_hessian",false);
  solver.setOption("tol",1e-18);
  solver.setOption("max_iter",1000);

  timer.start();
  solver.init();
  timer.stop("solver.init() ");

  VectorXd init_x;
  init_x.resize(G_Lame_rho.size());
  const vector<double> &g = tetmesh->material()._G;
  const vector<double> &la = tetmesh->material()._lambda;
  const vector<double> &rh = tetmesh->material()._rho;
  for (int i = 0; i < num_tet; ++i){
    init_x[i] = g[i];
    init_x[i+num_tet] = la[i];
    init_x[i+num_tet*2] = rh[i];
  }
  solver.setInput(&init_x[0],CasADi::NLP_X_INIT);

  cout << endl << endl;
  INFO_LOG("initial x: ");
  for (int i = 0; i < init_x.size(); ++i){
    cout << init_x[i] << " ,";
  }
  cout << endl << endl;

  vector<double> lower;
  lower.resize(G_Lame_rho.size());
  for (int i = 0; i < num_tet; ++i){
    lower[i] = 0.0f;
    lower[i+num_tet] = 0.0f;
    lower[i+num_tet*2] = 0.0f;
  }
  solver.setInput(lower,CasADi::NLP_LBX);

  // solving
  INFO_LOG("solving");
  timer.start();
  solver.solve();
  timer.stop("solver.solve() ");
  std::vector<double> vx(G_Lame_rho.size());
  solver.getOutput(vx,CasADi::NLP_X_OPT);  

  // reset the mtl for constrained elements
  for (int i = 0; i < num_tet*2; ++i){
    if (vx[i] > 1000.0f)
	  vx[i] = 1e-10;
  }

  // convert to Young's(E) and Poisson(v)
  vector<double> Young_E(num_tet), Poisson_v(num_tet);
  for (int i = 0; i < num_tet; ++i){

	const double gi = vx[i];
	const double li = vx[i+num_tet];
	const Matrix<double,2,1> Ev = ElasticMaterial<double>::fromLameConstant(gi,li);
	Young_E[i] = Ev(0,0);
	Poisson_v[i] = Ev(1,0);
	
	tetmesh->material()._G[i] = gi;
	tetmesh->material()._lambda[i] = li;
	tetmesh->material()._rho[i] = vx[i+num_tet*2];
  }

  // save results
  INFO_LOG("save results");

  cout << endl << endl << "G, L: \n";
  for (int i = 0; i < num_tet; ++i){
	cout << vx[i] << ", " << vx[i+num_tet] << ", " << vx[i+num_tet*2] << endl;
  }
  cout << endl << endl;

  cout << endl << endl << "E, v: \n";
  for (int i = 0; i < num_tet; ++i){
	cout << Young_E[i] << ", " << Poisson_v[i] << endl;
  }
  cout << endl << endl;

  bool succ = tetmesh->writeElasticMtlVTK(save_to);
  assert(succ);
}


void recoverFullMtlMethod3(pTetMesh tetmesh, const set<int> &fixednodes,
						   const MatrixXd &eigen_W, const VectorXd &eigen_lambda,
						   const double mu_stiff,const double mu_mass,const double mu_neigh,
   				           const string save_to){

  // load data.
  INFO_LOG("loading data");
  const int num_tet = (int)tetmesh->tets().size();
  const int n = tetmesh->nodes().size();

  SXMatrix W, lambda;
  CASADI::convert(eigen_W, W);
  assert_eq(W.size1(), n*3);
  assert_lt(W.size2(), W.size1());

  lambda.resize(eigen_lambda.size(), eigen_lambda.size());
  for (int i = 0; i < eigen_lambda.size(); ++i){
    lambda.elem(i,i) = eigen_lambda[i];
  }
  assert_eq(lambda.size1(), W.size2());
  assert_eq(lambda.size2(), W.size2());
  
  // init variables
  INFO_LOG("init variables");
  const CASADI::VSX G = CASADI::makeSymbolic(num_tet, "G");
  const CASADI::VSX Lame = CASADI::makeSymbolic(num_tet, "lame");
  const CASADI::VSX rho = CASADI::makeSymbolic(num_tet, "rho");
  const CASADI::VSX G_Lame_rho = CASADI::connect(CASADI::connect(G,Lame),rho);

  // compute K
  INFO_LOG("compute K");
  ComputeStiffnessMat elas(tetmesh);
  elas.setMaterial(G, Lame);
  assert(elas.prepare());
  VectorXd x0(n*3);
  tetmesh->nodes(x0);

  Timer timer;
  timer.start();
  SXMatrix K = elas.K(x0);
  K *= -1.0f;
  timer.stop("elas.K(x0) ");

  // compute mass matrix
  INFO_LOG("compute mass matrix");
  ComputeMassMat mass;
  SXMatrix M;
  mass.compute(M,*tetmesh,rho);

  /// remove fixed nodes
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

  // assemble objective function.
  INFO_LOG("assemble objective function.");

  SX objfun = 0.0f;
  const SXMatrix M1 = trans(W).mul(K.mul(W))-lambda;
  const SXMatrix M2 = trans(W).mul(M.mul(W))-CASADI::convert(MatrixXd((VectorXd::Ones(W.size2())).asDiagonal()));
  for (int i = 0; i < M1.size1(); ++i){
    for (int j = 0; j < M1.size2(); ++j){
	  objfun += mu_stiff*M1.elem(i,j)*M1.elem(i,j);
	  objfun += mu_mass*M2.elem(i,j)*M2.elem(i,j);
	}
  }

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
		objfun += 0.25f*mu_neigh*(vi+vj)*(G[i]-G[j])*(G[i]-G[j]);
		objfun += 0.25f*mu_neigh*(vi+vj)*(Lame[i]-Lame[j])*(Lame[i]-Lame[j]);
		objfun += 0.25f*mu_neigh*(vi+vj)*(rho[i]-rho[j])*(rho[i]-rho[j]);
	  }
	}
  }
  objfun *= 0.5f;

  // solve.
  INFO_LOG("init solver");
  CasADi::SXFunction fun = CasADi::SXFunction(G_Lame_rho,objfun);
  CasADi::IpoptSolver solver = CasADi::IpoptSolver(fun);

  solver.setOption("generate_hessian",false);
  solver.setOption("tol",1e-18);
  solver.setOption("max_iter",1000);

  timer.start();
  solver.init();
  timer.stop("solver.init() ");

  VectorXd init_x;
  init_x.resize(G_Lame_rho.size());
  const vector<double> &g = tetmesh->material()._G;
  const vector<double> &la = tetmesh->material()._lambda;
  const vector<double> &rh = tetmesh->material()._rho;
  for (int i = 0; i < num_tet; ++i){
    init_x[i] = g[i];
    init_x[i+num_tet] = la[i];
    init_x[i+num_tet*2] = rh[i];
  }
  solver.setInput(&init_x[0],CasADi::NLP_X_INIT);

  cout << endl << endl;
  INFO_LOG("initial x: ");
  for (int i = 0; i < init_x.size(); ++i){
    cout << init_x[i] << " ,";
  }
  cout << endl << endl;

  vector<double> lower;
  lower.resize(G_Lame_rho.size());
  for (int i = 0; i < num_tet; ++i){
    lower[i] = 0.0f;
    lower[i+num_tet] = 0.0f;
    lower[i+num_tet*2] = 0.0f;
  }
  solver.setInput(lower,CasADi::NLP_LBX);

  // solving
  INFO_LOG("solving");
  timer.start();
  solver.solve();
  timer.stop("solver.solve() ");
  std::vector<double> vx(G_Lame_rho.size());
  solver.getOutput(vx,CasADi::NLP_X_OPT);  

  // convert to Young's(E) and Poisson(v)
  vector<double> Young_E(num_tet), Poisson_v(num_tet);
  for (int i = 0; i < num_tet; ++i){

	const double gi = vx[i];
	const double li = vx[i+num_tet];
	const Matrix<double,2,1> Ev = ElasticMaterial<double>::fromLameConstant(gi,li);
	Young_E[i] = Ev(0,0);
	Poisson_v[i] = Ev(1,0);
	
	tetmesh->material()._G[i] = gi;
	tetmesh->material()._lambda[i] = li;
	tetmesh->material()._rho[i] = vx[i+num_tet*2];
  }

  // save results
  INFO_LOG("save results");

  cout << endl << endl << "G, L: \n";
  for (int i = 0; i < num_tet; ++i){
	cout << vx[i] << ", " << vx[i+num_tet] << ", " << vx[i+num_tet*2] << endl;
  }
  cout << endl << endl;

  cout << endl << endl << "E, v: \n";
  for (int i = 0; i < num_tet; ++i){
	cout << Young_E[i] << ", " << Poisson_v[i] << endl;
  }
  cout << endl << endl;

  bool succ = tetmesh->writeElasticMtlVTK(save_to);
  assert(succ);
}

void recoverCorrect(){

  // const string data_root = "/home/simba/Workspace/SolidSimulator/data/one_tet/model/";
  // const string data_root = "/home/simba/Workspace/SolidSimulator/data/two_tets/model/";
  // const string data_root = "/home/simba/Workspace/SolidSimulator/data/beam/model/";
  const string data_root = "/home/simba/Workspace/SolidSimulator/data/beam-coarse/model/";
  // const string data_root = "/home/simba/Workspace/SolidSimulator/data/beam-middle/model/";

  vector<int> fixednodes_vec;
  bool succ = loadVec(data_root+"/con_nodes.bou",fixednodes_vec, TEXT); assert(succ);
  set<int> fixednodes;
  for (int i = 0; i < fixednodes_vec.size(); ++i){
    fixednodes.insert(fixednodes_vec[i]);
  }

  computeEigenValues(data_root, 10, fixednodes);

  pTetMesh tetmesh = pTetMesh(new TetMesh());
  succ = tetmesh->load(data_root+"mesh.abq"); assert(succ);
  succ = tetmesh->loadElasticMtl(data_root+"tempt_mesh.elastic"); assert(succ);

  MatrixXd eigen_W;
  succ = load(data_root+"tempt_eigenvectors.b",eigen_W); assert(succ);

  VectorXd eigen_lambda;
  succ = load(data_root+"tempt_eigenvalues.b",eigen_lambda); assert(succ);

  //recoverFullMtl(tetmesh, fixednodes,eigen_W,eigen_lambda,1e-2,"./tempt/material_sim_3");
  // recoverFullMtl(tetmesh, fixednodes,eigen_W,eigen_lambda,1e-6,"./tempt/material_sim_2");
  // recoverFullMtl(tetmesh, fixednodes,eigen_W,eigen_lambda,1e-8,"./tempt/material_sim_1");
  // recoverFullMtl(tetmesh,fixednodes,eigen_W,eigen_lambda,0.000,"./tempt/material_sim_0");
  // recoverFullMtlMethod2(tetmesh,fixednodes,eigen_W,eigen_lambda,1e-6,"./tempt/material_sim_4");
  // recoverFullMtlMethod2(tetmesh,fixednodes,eigen_W,eigen_lambda,1e-2,"./tempt/material_sim_5");
  // recoverFullMtlMethod2(tetmesh,fixednodes,eigen_W,eigen_lambda,0.1,"./tempt/material_sim_6");
  // recoverFullMtlMethod3(tetmesh,fixednodes,eigen_W,eigen_lambda,0.1,"./tempt/material_sim_7");
  // recoverFullMtlMethod3(tetmesh,fixednodes,eigen_W,eigen_lambda,100,"./tempt/material_sim_8");
  // recoverFullMtlMethod3(tetmesh,fixednodes,eigen_W,eigen_lambda,0.0001,"./tempt/material_sim_9");
  // recoverFullMtlMethod3(tetmesh,fixednodes,eigen_W,eigen_lambda,1e-8,"./tempt/material_sim_10");
  // recoverFullMtlMethod3(tetmesh,fixednodes,eigen_W,eigen_lambda,1e-16,"./tempt/material_sim_16");
  recoverFullMtlMethod3(tetmesh,fixednodes,eigen_W,eigen_lambda,1.0,1.0,0.0f,"./tempt/material_sim_20");

}

void recoverOptMtl(){
  
  const string data_root = "/home/simba/Workspace/OthersProjects/AnimationEditor_sig/Data/beam_fine/model/";
  pTetMesh tetmesh = pTetMesh(new TetMesh());
  bool succ = tetmesh->load(data_root+"mesh.abq"); assert(succ);

  const string data_root2="/home/simba/Workspace/OthersProjects/AnimationEditor_sig/Data/beam_setmtl/model/";
  succ = tetmesh->loadElasticMtl(data_root2+"mesh.elastic"); assert(succ);
  // succ = tetmesh->loadElasticMtl(data_root+"mesh.elastic"); assert(succ);
  
  // tetmesh->writeElasticMtlVTK("./tempt/material_opt_ref");
  // exit(0);
  
  vector<int> fixednodes_vec;
  succ = loadVec(data_root+"/con_nodes.bou",fixednodes_vec, TEXT); assert(succ);
  set<int> fixednodes;
  for (int i = 0; i < fixednodes_vec.size(); ++i){
    fixednodes.insert(fixednodes_vec[i]);
  }

  MatrixXd eigen_W, S;
  succ = load(data_root+"scaled_W.b",eigen_W); assert(succ);
  succ = load(data_root+"S_opt_rlst.b",S); assert(succ);
  const MatrixXd t = eigen_W;
  eigen_W = t.leftCols(S.rows())*S;

  VectorXd eigen_lambda;
  succ = load(data_root+"lambda_opt_rlst.b",eigen_lambda); assert(succ);

  // succ = write("eigen_W.b",eigen_W); assert(succ);
  // succ = write("eigen_Lambda.b",eigen_lambda); assert(succ);
  // recoverFullMtl(tetmesh, fixednodes, eigen_W,eigen_lambda,1e12,"./tempt/material_opt_4");
  // recoverFullMtl(tetmesh, fixednodes, eigen_W,eigen_lambda,1e8,"./tempt/material_opt_5");
  // recoverFullMtl(tetmesh, fixednodes, eigen_W,eigen_lambda,10000,"./tempt/material_opt_3");
  // recoverFullMtl(tetmesh, fixednodes, eigen_W,eigen_lambda,100,"./tempt/material_opt_2");
  // recoverFullMtl(tetmesh, fixednodes, eigen_W,eigen_lambda,1.0,"./tempt/material_opt_1");
  // recoverFullMtl(tetmesh,fixednodes,eigen_W,eigen_lambda,0.001,"./tempt/material_opt_01");
  // recoverFullMtl(tetmesh,fixednodes,eigen_W,eigen_lambda,0.0,"./tempt/material_opt_0");
  // recoverFullMtl(tetmesh,fixednodes,eigen_W,eigen_lambda,1e-8,"./tempt/material_opt_1_8");
  // recoverFullMtl(tetmesh,fixednodes,eigen_W,eigen_lambda,1e-6,"./tempt/material_opt_1_6");

  // recoverFullMtlMethod2(tetmesh,fixednodes,eigen_W,eigen_lambda,1e-16,"./tempt/material_opt_24");
  // recoverFullMtlMethod2(tetmesh,fixednodes,eigen_W,eigen_lambda,1.0,"./tempt/material_opt_21");
  // recoverFullMtlMethod2(tetmesh,fixednodes,eigen_W,eigen_lambda,100,"./tempt/material_opt_22");
  // recoverFullMtlMethod2(tetmesh,fixednodes,eigen_W,eigen_lambda,0.001,"./tempt/material_opt_23");

  // recoverFullMtlMethod3(tetmesh,fixednodes,eigen_W,eigen_lambda,1e-16,"./tempt/material_opt_34");
  // recoverFullMtlMethod3(tetmesh,fixednodes,eigen_W,eigen_lambda,1.0,"./tempt/material_opt_31");
  // recoverFullMtlMethod3(tetmesh,fixednodes,eigen_W,eigen_lambda,100,"./tempt/material_opt_32");
  // recoverFullMtlMethod3(tetmesh,fixednodes,eigen_W,eigen_lambda,0.001,"./tempt/material_opt_33");
  // recoverFullMtlMethod3(tetmesh,fixednodes,eigen_W,eigen_lambda,10000,"./tempt/material_opt_35");

  // recoverFullMtlMethod3(tetmesh,fixednodes,eigen_W,eigen_lambda,0.0,1.0,1.0,"./tempt/material_opt_36");
  // recoverFullMtlMethod3(tetmesh,fixednodes,eigen_W,eigen_lambda,1.0,0.0,1.0,"./tempt/material_opt_37");
  // recoverFullMtlMethod3(tetmesh,fixednodes,eigen_W,eigen_lambda,1.0,1.0,1.0,"./tempt/material_opt_38");

  recoverFullMtlMethod3(tetmesh,fixednodes,eigen_W,eigen_lambda,1.0,1.0,10000.0,"./tempt/material_opt_40");
  // recoverFullMtlMethod3(tetmesh,fixednodes,eigen_W,eigen_lambda,1.0,1.0,1.0,"./tempt/material_opt_31");
  // recoverFullMtlMethod3(tetmesh,fixednodes,eigen_W,eigen_lambda,1.0,1.0,100.0,"./tempt/material_opt_32");

  // recoverFullMtlMethod2(tetmesh,fixednodes,eigen_W,eigen_lambda,0.0,"./tempt/material_opt_20");
  // recoverFullMtlMethod2(tetmesh,fixednodes,eigen_W,eigen_lambda,1.0,"./tempt/material_opt_21");
  // recoverFullMtlMethod2(tetmesh,fixednodes,eigen_W,eigen_lambda,100.0,"./tempt/material_opt_22");

  // recoverFullMtlMethod3(tetmesh,fixednodes,eigen_W,eigen_lambda,1.0,1.0,0.0,"./tempt/material_opt_50");
  // recoverFullMtlMethod3(tetmesh,fixednodes,eigen_W,eigen_lambda,1.0,1.0,1.0,"./tempt/material_opt_51");
  // recoverFullMtlMethod3(tetmesh,fixednodes,eigen_W,eigen_lambda,1.0,1.0,100.0,"./tempt/material_opt_52");
}

void recoverOptMtlCorrect(){
  
  const string data_root = "/home/simba/Workspace/AnimationEditor/Data/beam_fine/model/";
  pTetMesh tetmesh = pTetMesh(new TetMesh());
  bool succ = tetmesh->load(data_root+"mesh.abq"); assert(succ);

  succ = tetmesh->loadElasticMtl(data_root+"mesh.elastic"); assert(succ);
  
  vector<int> fixednodes_vec;
  succ = loadVec(data_root+"/con_nodes.bou",fixednodes_vec, TEXT); assert(succ);
  set<int> fixednodes;
  for (int i = 0; i < fixednodes_vec.size(); ++i){
    fixednodes.insert(fixednodes_vec[i]);
  }

  MatrixXd eigen_W;
  succ = load(data_root+"eigenvectors.b",eigen_W); assert(succ);
  const MatrixXd t = eigen_W;
  eigen_W = t.leftCols(2);

  VectorXd eigen_lambda;
  succ = load(data_root+"eigenvalues.b",eigen_lambda); assert(succ);
  const VectorXd te = eigen_lambda;
  eigen_lambda = te.head(2);

  recoverFullMtl(tetmesh, fixednodes, eigen_W, eigen_lambda, 1e6, "./tempt/material_opt_c");
}

int main(int argc, char *argv[]){

  // testK();
  // testM();
  // testKSMall();
  // recoverCorrect();
  recoverOptMtl();
  // recoverOptMtlCorrect();
}
