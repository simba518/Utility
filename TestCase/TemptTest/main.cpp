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
#include "ComputeStiffnessMat.h"
#include "ComputeMassMat.h"
#include "MaterialFitting.h"
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

void test_diagM(){
  
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

void testM(){
  
  const string tet_fname = std::string(TEST_DATA_DIR)+"two_tet.abq";
  const string tet_mtl_fname = std::string(TEST_DATA_DIR)+"two_tet.elastic";
  pTetMesh mesh = pTetMesh(new TetMesh);
  assert ( mesh->load(tet_fname) );
  assert ( mesh->loadElasticMtl(tet_mtl_fname) );
  const int n = mesh->nodes().size();

  MassMatrix mass_c;
  SparseMatrix<double> Me;
  mass_c.compute(Me, *mesh);

  ComputeMassMat mass;
  SXMatrix M;
  mass.compute(M, *mesh, false);
  MatrixXd Ms;
  CASADI::convert(M,Ms);

  assert_eq(Me.rows(), Me.cols());
  assert_eq(Me.rows(), Ms.rows());
  assert_eq(Me.rows(), Ms.cols());
  const MatrixXd mm = Ms-MatrixXd(Me);
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
  SparseMatrix<double> M;
  mass.compute(M,*tetmesh);

  // remove fixed nodes
  SparseMatrix<double> P;
  EIGEN3EXT::genReshapeMatrix(K.rows(),3,fixednodes,P);
  K = P*(K*P.transpose());
  const SparseMatrix<double> Msub = P*(M*P.transpose());
  const SparseMatrix<double> Mlower = EIGEN3EXT::getLower(Msub);

  /// compute W, lambda
  INFO_LOG("compute W, lambda");
  MatrixXd W_full;
  VectorXd lambda_full;
  const SparseMatrix<double> Klower = EIGEN3EXT::getLower(K);
  succ = EigenSparseGenEigenSolver::solve(Klower,Mlower,W_full,lambda_full,eigenNum); 
  assert(succ);

  // const MatrixXd W = W_full.rightCols(lambda_full.size()-6);
  // const VectorXd lambda = lambda_full.tail(lambda_full.size()-6);

  // test Wt*M*W
  INFO_LOG("check results");
  const MatrixXd WtMWI = W_full.transpose()*Msub*W_full-MatrixXd::Identity(W_full.cols(),W_full.cols());
  const MatrixXd KW_MWLambda = K*W_full-(Msub*W_full)*lambda_full.asDiagonal();
  cout << "(WtMW-I).norm(): " << WtMWI.norm() << endl;
  cout << "(KW_MWLambda).norm(): " << KW_MWLambda.norm() << endl;
  cout << "eigenvalues: " << lambda_full.transpose() << endl;
  cout << "norm(Klower): " << Klower.norm() << endl;
  cout << "norm(M): " << Msub.norm() << endl;

  // add fixed nodes to W, and lambda, then save
  const MatrixXd W = P.transpose()*W_full;
  const VectorXd lambda = lambda_full;
  succ = write(data_root+"tempt_eigenvalues.b", lambda); assert(succ);
  succ = write(data_root+"tempt_eigenvectors.b", W); assert(succ);
  
}

void recoverOpt(){
  
  const string data_root = "/home/simba/Workspace/AnimationEditor/Data/beam-coarse/model/";
  MatrixXd eig_W, S;
  bool succ = load(data_root+"scaled_W.b",eig_W); assert(succ);
  succ = load(data_root+"S.b",S); assert(succ);
  // S = MatrixXd::Identity(S.rows(),S.cols());
  const MatrixXd t = eig_W;
  eig_W = t.leftCols(S.rows())*S;
  
  VectorXd eig_lambda;
  succ = load(data_root+"lambda.b",eig_lambda); assert(succ);
  const double h = 0.03f;
  eig_lambda *= 1.0f/(h*h);

  succ = write("eig_lambda.b",eig_lambda); assert(succ);
  succ = write("eig_W.b",eig_W); assert(succ);

  const int used_r = 3;
  MatrixXd W = eig_W.leftCols(used_r);
  VectorXd lambda = eig_lambda.head(used_r);
  
  // MatrixXd W = eig_W.col(1);
  // VectorXd lambda = eig_lambda.segment<1>(1);

  cout<< "lambda: " << eig_lambda.transpose() << endl;

  MaterialFitting_Diag_M mtlfit_m;
  { // fit density
  	INFO_LOG("fit density");
  	mtlfit_m.loadTetMesh(data_root+"mesh.abq");
  	mtlfit_m.loadMtl(data_root+"mesh.elastic");
  	mtlfit_m.loadFixednodes(data_root+"/con_nodes.bou");
  	mtlfit_m.setWLambda(W, lambda);
  	mtlfit_m.computeK();
  	mtlfit_m.computeM();
  	mtlfit_m.removeFixedDOFs();
  	mtlfit_m.useHessian(true);

	mtlfit_m.setBounds(1e-20,1e8);
  	mtlfit_m.setMuSmoothGL(0, 0);
  	mtlfit_m.setMuSmoothEv(0, 0);
	mtlfit_m.setMuSmoothDensity(1.0);
  	mtlfit_m.setMuAverageDensity(0.0f);

  	mtlfit_m.assembleObjfun();
  	// mtlfit_m.solveByIpopt();
  	mtlfit_m.solveByLinearSolver();
  	mtlfit_m.saveResults("./tempt/material_opt_m_30");
	mtlfit_m.printResult();
  }
  
  { // transform W
  	INFO_LOG("transform W");
  	DiagonalMatrix<double,-1> inv_sqrt_M_real, sqrt_M_new;
  	mtlfit_m.computeM(inv_sqrt_M_real);
  	mtlfit_m.computeM(sqrt_M_new, mtlfit_m.getDensityResult());
  	const DiagonalMatrix<double,-1> M_real= inv_sqrt_M_real;
  	const DiagonalMatrix<double,-1> M_new = sqrt_M_new;
  	cout<< "\n\nW^t*M*W: \n" << (W.transpose()*M_new*W) << endl;
  }

  pMaterialFitting mtlfit_k = pMaterialFitting(new MaterialFitting_EV_MA_K);
  { // fit G, l
  	INFO_LOG("fit G, l");
  	mtlfit_k->loadTetMesh(data_root+"mesh.abq");
  	mtlfit_k->loadMtl(data_root+"mesh.elastic");
  	mtlfit_k->loadFixednodes(data_root+"/con_nodes.bou");
  	mtlfit_k->setWLambda(W, lambda);
  	mtlfit_k->computeK();
  	mtlfit_k->computeM();
  	mtlfit_k->removeFixedDOFs();
  	mtlfit_k->useHessian(true);

  	mtlfit_k->setBounds(1e-8,1000);
  	mtlfit_k->setMuSmoothGL(0, 0);
  	mtlfit_k->setMuSmoothEv(2e3, 0);
  	mtlfit_k->setMuSmoothDensity(0.0f);
  	mtlfit_k->setMuAverageDensity(0.0f);

  	mtlfit_k->assembleObjfun();
  	mtlfit_k->solveByIpopt();
  	// mtlfit_k->solveByLinearSolver();
  	mtlfit_k->saveResults("./tempt/material_opt_k_30");
  	mtlfit_k->printResult();
  }

  { // check resutls
	SparseMatrix<double> K, K0;
  	mtlfit_k->computeK(K,mtlfit_k->getShearGResult(),mtlfit_k->getLameResult());
  	mtlfit_k->computeK(K0);

	DiagonalMatrix<double,-1> M0_diag;
  	mtlfit_k->computeM(M0_diag);

	const SparseMatrix<double> P = mtlfit_k->getMatrixForRemovingFixedDOFs();
	K = P*K*P.transpose();
	K0 = P*K0*P.transpose();
	W = P*W;
	const SparseMatrix<double> M0 = P*M0_diag*P.transpose();

  	const MatrixXd la = lambda.asDiagonal();
  	const MatrixXd m0 = (W.transpose()*K0*W-la);
  	const MatrixXd m1 = (W.transpose()*K*W-la);
	const MatrixXd m2 = (W.transpose()*K0*W-(W.transpose()*M0*W)*la);
	const MatrixXd m3 = (W.transpose()*K*W-(W.transpose()*M0*W)*la);

  	cout << "\n\nlambda: " << lambda.transpose() << "\n\n";
  	cout << "\n\nW^t*K0*W-Lambda: "<<m0.norm()<<"\n\n"<<m0<<"\n\n";
  	cout << "\n\nW^t*K*W-Lambda:  "<<m1.norm()<<"\n\n"<<m1<<"\n\n";
  	cout << "\n\nW^t*K0*W-W^t*M0*W*Lambda:  "<<m2.norm()<<"\n\n"<<m2<<"\n\n";
  	cout << "\n\nW^t*K*W-W^t*M0*W*Lambda:  "<<m3.norm()<<"\n\n"<<m3<<"\n\n";
  }
}

void recoverSim(){

  const string data_root = "/home/simba/Workspace/SolidSimulator/data/beam-coarse/model/";
  const string fixed_nodes = data_root + "/con_nodes.bou";
  vector<int> fixednodes_vec;
  bool succ = loadVec(fixed_nodes,fixednodes_vec, TEXT); assert(succ);
  set<int> fixednodes;
  for (int i = 0; i < fixednodes_vec.size(); ++i)
    fixednodes.insert(fixednodes_vec[i]);
  computeEigenValues(data_root, 5, fixednodes);

  MatrixXd eigen_W;
  succ = load(data_root+"tempt_eigenvectors.b",eigen_W); assert(succ);
  VectorXd eigen_lambda;
  succ = load(data_root+"tempt_eigenvalues.b",eigen_lambda); assert(succ);
  cout<< "\n\nLambda: " << eigen_lambda.transpose() << "\n\n";
  cout << eigen_W.col(0).norm() << endl;
  // eigen_W.col(0) += VectorXd::Random(eigen_W.rows())*eigen_W.col(0).norm()*0.001f;
  // cout << eigen_W.col(0).norm() << endl;

  MaterialFitting_Diag_M mtlfit_m;
  { // fit density
	INFO_LOG("fit density");
	mtlfit_m.loadTetMesh(data_root+"mesh.abq");
	mtlfit_m.loadMtl(data_root+"tempt_mesh.elastic");
	mtlfit_m.loadFixednodes(fixed_nodes);
	mtlfit_m.setWLambda(eigen_W, eigen_lambda);
	mtlfit_m.computeK();
	mtlfit_m.computeM();
	mtlfit_m.removeFixedDOFs();
	mtlfit_m.useHessian(true);

  	mtlfit_m.setMuSmoothGL(0, 0);
  	mtlfit_m.setMuSmoothEv(0, 0);
	mtlfit_m.setMuSmoothDensity(1e-10);
  	mtlfit_m.setMuAverageDensity(0.0f);

	mtlfit_m.assembleObjfun();
	mtlfit_m.solveByIpopt();
	mtlfit_m.saveResults("./tempt/material_sim_m_30");
  }

  { // check H
	MatrixXd H;
	VectorXd g;
	mtlfit_m.hessGrad(H,g);
	SelfAdjointEigenSolver<MatrixXd> es;
	es.compute(H);
	cout << "H.norm() = " << H.norm() << endl;
	cout << "The eigenvalues of H are: \n" << es.eigenvalues().transpose() << endl;
  }

  { // check W^t*M*W
  	SparseMatrix<double> M;
  	mtlfit_m.computeM(M, mtlfit_m.getDensityResult());
	const int r = eigen_W.cols();
	cout <<"\n\nnorm:"<<(eigen_W.transpose()*M*eigen_W-MatrixXd::Identity(r,r)).norm()<<endl;
  	cout << "\n\nW^t*M*W: \n\n\n" << (eigen_W.transpose()*M*eigen_W) << endl;
  }
  
  // { // transform W
  // 	INFO_LOG("transform W");
  // 	SparseMatrix<double> M_real, M_new;
  // 	mtlfit_m.computeM(M_real);
  // 	mtlfit_m.computeM(M_new, mtlfit_m.getDensityResult());
	
  // 	const MatrixXd dM_real = M_real;
  // 	LLT<MatrixXd> lltOfM_real(dM_real);
  // 	const MatrixXd L_real = lltOfM_real.matrixL();

  // 	const MatrixXd dM_new = M_new;
  // 	LLT<MatrixXd> lltOfM_new(dM_new);
  // 	const MatrixXd L_new = lltOfM_new.matrixL();

  // 	eigen_W = L_real.transpose().lu().solve(L_new.transpose()*eigen_W);
  // 	cout<< "\n\nW^t*M*W: \n\n\n" << (eigen_W.transpose()*M_real*eigen_W) << endl;

  // 	// INFO_LOG("transform W");
  // 	// DiagonalMatrix<double,-1> inv_sqrt_M_real, sqrt_M_new;
  // 	// mtlfit_m.computeM(inv_sqrt_M_real);
  // 	// mtlfit_m.computeM(sqrt_M_new, mtlfit_m.getDensityResult());
  // 	// const DiagonalMatrix<double,-1> M_real= inv_sqrt_M_real;
  // 	// const DiagonalMatrix<double,-1> M_new = sqrt_M_new;
  // 	// cout<< "\n\nW^t*M*W: \n" << (eigen_W.transpose()*M_new*eigen_W) << endl;
  // 	// for (int i = 0; i < inv_sqrt_M_real.rows(); ++i){
  // 	//   assert_gt(inv_sqrt_M_real.diagonal()[i], 0);
  // 	//   assert_gt(sqrt_M_new.diagonal()[i], 0);
  // 	//   inv_sqrt_M_real.diagonal()[i] = 1.0f/sqrt(inv_sqrt_M_real.diagonal()[i]);
  // 	//   sqrt_M_new.diagonal()[i] = sqrt(sqrt_M_new.diagonal()[i]);
  // 	// }
  // 	// eigen_W = inv_sqrt_M_real*(sqrt_M_new*eigen_W);
  // 	// cout<< "\n\nW^t*M0*W: \n\n\n" << (eigen_W.transpose()*M_real*eigen_W) << endl;
  // }

  // MaterialFitting_EV_Diag_K mtlfit_k;
  // { // fit G, l
  // 	INFO_LOG("fit G, l");
  // 	mtlfit_k.loadTetMesh(data_root+"mesh.abq");
  // 	mtlfit_k.loadMtl(data_root+"tempt_mesh.elastic");
  // 	mtlfit_k.loadFixednodes(fixed_nodes);
  // 	mtlfit_k.setWLambda(eigen_W, eigen_lambda);
  // 	mtlfit_k.computeK();
  // 	mtlfit_k.computeM();
  // 	mtlfit_k.removeFixedDOFs();
  // 	mtlfit_k.useHessian(true);

  // 	mtlfit_k.setBounds(1e3,2e7);
  // 	mtlfit_k.setMuSmoothGL(0, 0);
  // 	mtlfit_k.setMuSmoothEv(1e-12, 0);
  // 	mtlfit_k.setMuSmoothDensity(0.0f);
  // 	mtlfit_k.setMuAverageDensity(0.0f);

  // 	mtlfit_k.assembleObjfun();
  // 	mtlfit_k.solveByIpopt();
  // 	// mtlfit_k.solveByLinearSolver();
  // 	// mtlfit_k.solveByNNLS();
  // 	mtlfit_k.saveResults("./tempt/material_sim_k_30");

  // 	SparseMatrix<double> K, K0;
  // 	mtlfit_k.computeK(K,mtlfit_k.getShearGResult(),mtlfit_k.getLameResult());
  // 	mtlfit_k.computeK(K0);
  // 	const MatrixXd la = eigen_lambda.asDiagonal();
  // 	const MatrixXd m1 = (eigen_W.transpose()*K*eigen_W-la);
  // 	const MatrixXd m0 = (eigen_W.transpose()*K0*eigen_W-la);
  // 	cout << "\n\nW^t*K0*W-Lambda:\n\n"<<m0.norm()<<"\n\n"<<m0<<"\n\n";
  // 	cout << "\n\nW^t*K*W-Lambda:\n\n"<<m1.norm()<<"\n\n"<<m1<<"\n\n";
  // 	mtlfit_k.printResult();
  // }
}

void recoverSim_MA_K(){

  const string data_root = "/home/simba/Workspace/SolidSimulator/data/beam-coarse/model/";
  const string fixed_nodes = data_root + "/con_nodes.bou";
  vector<int> fixednodes_vec;
  bool succ = loadVec(fixed_nodes,fixednodes_vec, TEXT); assert(succ);
  set<int> fixednodes;
  for (int i = 0; i < fixednodes_vec.size(); ++i)
    fixednodes.insert(fixednodes_vec[i]);
  computeEigenValues(data_root, 4, fixednodes);

  MatrixXd eigen_W;
  succ = load(data_root+"tempt_eigenvectors.b",eigen_W); assert(succ);
  VectorXd eigen_lambda;
  succ = load(data_root+"tempt_eigenvalues.b",eigen_lambda); assert(succ);

  MaterialFitting_MA_K mtlfit;
  { // fit density
	INFO_LOG("fit density");
	mtlfit.loadTetMesh(data_root+"mesh.abq");
	mtlfit.loadMtl(data_root+"tempt_mesh.elastic");
	mtlfit.loadFixednodes(fixed_nodes);
	mtlfit.setWLambda(eigen_W, eigen_lambda);
	mtlfit.computeK();
	mtlfit.computeM();
	mtlfit.removeFixedDOFs();
	mtlfit.useHessian(true);

  	mtlfit.setMuSmoothGL(1e-15, 1e-15);
  	mtlfit.setMuSmoothEv(0, 0);
	mtlfit.setMuSmoothDensity(0.0f);
  	mtlfit.setMuAverageDensity(0.0f);

	mtlfit.assembleObjfun();
	mtlfit.solveByIpopt();
	// mtlfit.solveByNNLS();
	mtlfit.saveResults("./tempt/material_sim_30");
  }
}

void recoverOpt_MA_K(){

  const string data_root = "/home/simba/Workspace/AnimationEditor/Data/beam-coarse/model/";
  MatrixXd eigen_W, S;
  bool succ = load(data_root+"scaled_W.b",eigen_W); assert(succ);
  succ = load(data_root+"S_opt_rlst.b",S); assert(succ);
  const MatrixXd t = eigen_W;
  eigen_W = t.leftCols(S.rows())*S;

  VectorXd eigen_lambda;
  succ = load(data_root+"lambda_opt_rlst.b",eigen_lambda); assert(succ);

  MaterialFitting_MA_K mtlfit;
  { // fit density
  	INFO_LOG("fit density");
  	mtlfit.loadTetMesh(data_root+"mesh.abq");
  	mtlfit.loadMtl(data_root+"mesh.elastic");
  	mtlfit.loadFixednodes(data_root+"/con_nodes.bou");
  	mtlfit.setWLambda(eigen_W.leftCols(4), eigen_lambda.head(4));
  	mtlfit.computeK();
  	mtlfit.computeM();
  	mtlfit.removeFixedDOFs();
  	mtlfit.useHessian(true);

  	mtlfit.setMuSmoothGL(1000, 1000);
  	mtlfit.setMuSmoothEv(0, 0);
	mtlfit.setMuSmoothDensity(0.0f);
  	mtlfit.setMuAverageDensity(0.0f);
  	mtlfit.setMuAverageDensity(0.0f);

  	mtlfit.assembleObjfun();
  	mtlfit.solveByIpopt();
  	// mtlfit.solveByNNLS();
  	mtlfit.saveResults("./tempt/material_opt_m_30");
  }
}

int main(int argc, char *argv[]){

  // testK();
  // testKSMall();
  // testM();
  // test_diagM();
  // recoverSim();
  recoverOpt();
  // recoverSim_MA_K();
  // recoverOpt_MA_K();
}

// { // transform W
// 	INFO_LOG("transform W");
// 	SparseMatrix<double> M_real, M_new;
// 	mtlfit_m.computeM(M_real);
// 	mtlfit_m.computeM(M_new, mtlfit_m.getDensityResult());
	
// 	const MatrixXd dM_real = M_real;
// 	LLT<MatrixXd> lltOfM_real(dM_real);
// 	const MatrixXd L_real = lltOfM_real.matrixL();

// 	const MatrixXd dM_new = M_new;
// 	LLT<MatrixXd> lltOfM_new(dM_new);
// 	const MatrixXd L_new = lltOfM_new.matrixL();

// 	const MatrixXd W_old = W;
// 	W = L_real.transpose().lu().solve(L_new.transpose()*W);
// 	cout<< "\n\nW^t*M*W: \n\n\n" << (W.transpose()*M_real*W) << endl;
// 	cout<< "W_old: " << W_old.norm() << "\n";
// 	cout<< "W_new: " << W.norm() << "\n\n";
// }
  


// for (int i = 0; i < inv_sqrt_M_real.rows(); ++i){
//   assert_gt(inv_sqrt_M_real.diagonal()[i], 0);
//   assert_gt(sqrt_M_new.diagonal()[i], 0);
//   inv_sqrt_M_real.diagonal()[i] = 1.0f/sqrt(inv_sqrt_M_real.diagonal()[i]);
//   sqrt_M_new.diagonal()[i] = sqrt(sqrt_M_new.diagonal()[i]);
// }
// const MatrixXd W1 = inv_sqrt_M_real*(sqrt_M_new*W);
// const VectorXd l1 = lambda;
// W = W1.col(0);
// lambda = l1.head(1);
// cout<< "\n\nW^t*M0*W: \n\n\n" << (W.transpose()*M_real*W) << endl;
