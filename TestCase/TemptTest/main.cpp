#include <iostream>
#include <MatrixTools.h>
#include <MatrixIO.h>
#include <interfaces/ipopt/ipopt_solver.hpp>
#include <CASADITools.h>
#include <MassMatrix.h>
#include <SparseGenEigenSolver.h>
#include <SparseMatrixTools.h>
#include <ElasticForceTetFullStVK.h>
#include "ComputeStiffnessMat.h"
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
  elas.K(x0);
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

void computeEigenValues(){
  
  const string data_root = "/home/simba/Workspace/SolidSimulator/data/one_tet/model/";
  
  // load data.
  pTetMesh tetmesh = pTetMesh(new TetMesh());
  bool succ = tetmesh->load(data_root+"mesh.abq"); assert(succ);
  succ = tetmesh->loadElasticMtl(data_root+"mesh.elastic"); assert(succ);
  const int num_tet = (int)tetmesh->tets().size();
  const int n = tetmesh->nodes().size();

  // compute K
  ElasticForceTetFullStVK elas(tetmesh);
  assert(elas.prepare());
  VectorXd x0(n*3);
  tetmesh->nodes(x0);
  SparseMatrix<double> K = elas.K(x0);
  K *= -1.0f;

  // compute mass matrix
  MassMatrix mass;
  DiagonalMatrix<double,-1> diagM;
  mass.compute(diagM,*tetmesh);

  /// compute W, lambda
  MatrixXd W;
  VectorXd lambda;
  const int eigenNum = 11;
  const SparseMatrix<double> Klower = EIGEN3EXT::getLower(K);
  succ = EigenSparseGenEigenSolver::solve(Klower,diagM,W,lambda,eigenNum); assert(succ);

  // test Wt*M*W
  const MatrixXd WtMWI = W.transpose()*diagM*W-MatrixXd::Identity(W.cols(),W.cols());
  const MatrixXd KW_MWLambda = K*W-(diagM*W)*lambda.asDiagonal();
  cout << "(WtMW-I).norm(): " << WtMWI.norm() << endl;
  cout << "(KW_MWLambda).norm(): " << KW_MWLambda.norm() << endl;
  cout << "eigenvalues: " << lambda.transpose() << endl;
  cout << "norm(Klower): " << Klower.norm() << endl;
  cout << "norm(M): " << diagM.diagonal().norm() << endl;

  succ = write(data_root+"tempt_eigenvalues.b", lambda); assert(succ);
  succ = write(data_root+"tempt_eigenvectors.b", W); assert(succ);
  
}

void recoverFullMtl(){

  const string data_root = "/home/simba/Workspace/SolidSimulator/data/one_tet/model/";
  
  // load data.
  INFO_LOG("loading data");
  pTetMesh tetmesh = pTetMesh(new TetMesh());
  bool succ = tetmesh->load(data_root+"mesh.abq"); assert(succ);
  succ = tetmesh->loadElasticMtl(data_root+"tempt_mesh.elastic"); assert(succ);
  const int num_tet = (int)tetmesh->tets().size();
  const int n = tetmesh->nodes().size();

  SXMatrix W, lambda;
  MatrixXd eigen_W;
  succ = load(data_root+"tempt_eigenvectors.b",eigen_W); assert(succ);
  CASADI::convert(eigen_W, W);
  assert_eq(W.size1(), n*3);
  assert_lt(W.size2(), W.size1());

  VectorXd eigen_lambda;
  succ = load(data_root+"tempt_eigenvalues.b",eigen_lambda); assert(succ);
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
  elas.K(x0);
  SXMatrix K = elas.K(x0);
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

  // assemble objective function.
  INFO_LOG("assemble objective function.");
  const double mu_G = 0.0f;
  const double mu_L = 0.0f;

  SX objfun = 0.0f;
  const SXMatrix M1 = K.mul(W)-M.mul(W.mul(lambda));
  for (int i = 0; i < M1.size1(); ++i){
    for (int j = 0; j < M1.size2(); ++j)
	 objfun += M1.elem(i,j)*M1.elem(i,j);
  }

  const vector<double> &g = tetmesh->material()._G;
  const vector<double> &la = tetmesh->material()._lambda;
  assert_eq(g.size(), num_tet);
  assert_eq(la.size(), num_tet);
  for (int i = 0; i < num_tet; ++i){
	objfun += mu_G*(G[i]-g[i])*(G[i]-g[i]);
	objfun += mu_L*(Lame[i]-la[i])*(Lame[i]-la[i]);
  }
  objfun *= 0.5f;

  // solve.
  INFO_LOG("init solver");
  CasADi::SXFunction fun = CasADi::SXFunction(G_Lame,objfun);
  CasADi::IpoptSolver solver = CasADi::IpoptSolver(fun);

  solver.setOption("generate_hessian",true);
  solver.setOption("tol",1e-8);
  solver.setOption("max_iter",100);
  solver.init();

  VectorXd init_x;
  init_x.resize(G_Lame.size());
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

  INFO_LOG("solving");
  solver.solve();

  // save results.
  INFO_LOG("save results");
  std::vector<double> vx(G_Lame.size());
  solver.getOutput(vx,CasADi::NLP_X_OPT);  

  cout << endl << endl;
  for (int i = 0; i < vx.size(); ++i){
    cout << vx[i] << " ,";
  }
  cout << endl << endl;
  
}

int main(int argc, char *argv[]){

  // testK();
  computeEigenValues();
  recoverFullMtl();
}
