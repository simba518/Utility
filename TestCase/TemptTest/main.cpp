#include <MatrixTools.h>
#include <MatrixIO.h>
#include <SparseGenEigenSolver.h>
#include <SparseMatrixTools.h>
#include <Timer.h>
#include "MaterialFitting.h"
using namespace UTILITY;
using namespace EIGEN3EXT;
using namespace ELASTIC_OPT;

void computeEigenValues(const string data_root, const string mtl_file,
						const int eigenNum,const set<int> &fixednodes){
  
  // load data.
  INFO_LOG("load data");
  pTetMesh tetmesh = pTetMesh(new TetMesh());
  bool succ = tetmesh->load(data_root+"mesh.abq"); assert(succ);
  succ = tetmesh->loadElasticMtl(mtl_file); assert(succ);
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

void recover(const string data_root, const string mtl_file, const string save_to,
			 const MatrixXd &W, const VectorXd &lambda, 
			 const double Gs=1e-3, const double Gl=1e-3,
			 const double lower=1, const double upper=1e10){

  MaterialFitting_MA_K mtlfit_k;

  { // fit G, l
  	INFO_LOG("fit G, l");
  	mtlfit_k.loadTetMesh(data_root+"mesh.abq");
  	mtlfit_k.loadMtl(mtl_file);
  	mtlfit_k.loadFixednodes(data_root+"/con_nodes.bou");
  	mtlfit_k.setWLambda(W, lambda);
  	mtlfit_k.computeK();
  	mtlfit_k.computeM();
  	mtlfit_k.removeFixedDOFs();
  	mtlfit_k.useHessian(false);

  	mtlfit_k.setBounds(lower,upper);
  	mtlfit_k.setMuSmoothGL(Gs, Gl);
  	mtlfit_k.setMuSmoothEv(0, 0);
  	mtlfit_k.setMuSmoothDensity(0.0f);
  	mtlfit_k.setMuAverageDensity(0.0f);
  	mtlfit_k.assembleObjfun();

  	// mtlfit_k.solveByIpopt();
	// mtlfit_k.solveByAlglib();
	mtlfit_k.solveByMPRGP();

  	mtlfit_k.saveResults(save_to);
  	mtlfit_k.printResult();
  }

  { // check resutls
	SparseMatrix<double> K, K0;
  	mtlfit_k.computeK(K,mtlfit_k.getShearGResult(),mtlfit_k.getLameResult());
  	mtlfit_k.computeK(K0);

	DiagonalMatrix<double,-1> M0_diag;
  	mtlfit_k.computeM(M0_diag);

	const SparseMatrix<double> P = mtlfit_k.getMatrixForRemovingFixedDOFs();
	K = P*K*P.transpose();
	K0 = P*K0*P.transpose();
	const MatrixXd W_sub = P*W;
	const SparseMatrix<double> M0 = P*M0_diag*P.transpose();

  	const MatrixXd la = lambda.asDiagonal();
  	const MatrixXd m0 = (W_sub.transpose()*K0*W_sub-la);
  	const MatrixXd m1 = (W_sub.transpose()*K*W_sub-la);
	const MatrixXd m2 = (W_sub.transpose()*K0*W_sub-(W_sub.transpose()*M0*W_sub)*la);
	const MatrixXd m3 = (W_sub.transpose()*K*W_sub-(W_sub.transpose()*M0*W_sub)*la);

  	cout << "\n\nlambda: " << lambda.transpose() << "\n\n";
  	cout << "\n\nW^t*K0*W-Lambda: "<<m0.norm()<<"\n\n"<<m0<<"\n\n";
  	cout << "\n\nW^t*K*W-Lambda:  "<<m1.norm()<<"\n\n"<<m1<<"\n\n";
  	cout << "\n\nW^t*K0*W-W^t*M0*W*Lambda:  "<<m2.norm()<<"\n\n"<<m2<<"\n\n";
  	cout << "\n\nW^t*K*W-W^t*M0*W*Lambda:  "<<m3.norm()<<"\n\n"<<m3<<"\n\n";
	cout << "W^t*M0*W: \n\n" << W_sub.transpose()*M0*W_sub << "\n\n";
  }

}

void recoverSim(const string data_root, const string data_root_opt, const int used_r, const double h){

  { // compute eigenvalue
	const string fixed_nodes = data_root + "/con_nodes.bou";
	vector<int> fixednodes_vec;
	bool succ = loadVec(fixed_nodes,fixednodes_vec, TEXT); assert(succ);
	set<int> fixednodes;
	for (int i = 0; i < fixednodes_vec.size(); ++i)
	  fixednodes.insert(fixednodes_vec[i]);
	computeEigenValues(data_root, data_root+"mesh.elastic", 5, fixednodes);
  }

  MatrixXd eig_W;
  VectorXd eig_lambda;

  { // load opt value
	// load W
	MatrixXd S;
	bool succ = load(data_root_opt+"scaled_W.b",eig_W); assert(succ);
	succ = load(data_root_opt+"S.b",S); assert(succ);
	const MatrixXd t = eig_W;
	eig_W = t.leftCols(S.rows())*S;
	
	// load lambda
	succ = load(data_root_opt+"lambda.b",eig_lambda); assert(succ);
	eig_lambda *= 1.0f/(h*h);
  }

  { // load real value
	// bool succ = load(data_root+"tempt_eigenvectors.b",eig_W); assert(succ);
	// succ = load(data_root+"tempt_eigenvalues.b",eig_lambda); assert(succ);
  }

  const MatrixXd W = eig_W.leftCols(used_r);
  const VectorXd lambda = eig_lambda.head(used_r);
  cout<< "\n\nLambda: " << lambda.transpose() << "\n\n";
  recover(data_root, data_root_opt+"/mesh.elastic", "./tempt/test", W, lambda);

}

void recoverCoarse(){

  const string data_root = "/home/simba/Workspace/SolidSimulator/data/beam-coarse/model/";
  const string data_root_opt ="/home/simba/Workspace/AnimationEditor/Data/beam-coarse/model/";
  recoverSim(data_root, data_root_opt, 2, 0.03);
}

void recoverFine(){

  const string data_root = "/home/simba/Workspace/SolidSimulator/data/beam/";

  MatrixXd eig_W;
  VectorXd eig_lambda;

  bool succ = true;
  { // load opt value
	succ = load(data_root+"tempt/opt_W.b",eig_W); assert(succ);
	succ = load(data_root+"tempt/opt_lambda.b",eig_lambda); assert(succ);
  }

  { // load real value
	// succ = load(data_root+"tempt/eigenvectors_smooth.b",eig_W); assert(succ);
	// succ = load(data_root+"tempt/eigenvalues_smooth.b",eig_lambda); assert(succ);
  }

  MatrixXd W(eig_W.rows(),2);
  W.col(0) = eig_W.col(0);
  W.col(1) = eig_W.col(1);
  VectorXd lambda(2);
  lambda << eig_lambda[0], eig_lambda[1];
  cout<< "\n\nLambda: " << lambda.transpose() << "\n\n";

  /// G,l \approx E/3
  // solve
  recover(data_root+"model/",data_root+"model/mesh.elastic","./tempt/test2",
		  W,lambda,1e-8,1e-8,1e-3,1e16);

  recover(data_root+"model/",data_root+"model/mesh.elastic","./tempt/test1",
		  W,lambda,1e-8,1e-8,1e3,3e6);

  recover(data_root+"model/",data_root+"model/mesh.elastic","./tempt/test3",
		  W,lambda,1e-3,1e-3,1e-3,1e16);

  recover(data_root+"model/",data_root+"model/tempt_mesh.elastic","./tempt/test4",
		  W,lambda,1e-3,1e-3,1e-3,1e16);

  W.resize(W.rows(), 1);
  W.col(0) = eig_W.col(0);
  lambda.resize(1);
  lambda << eig_lambda[0];
  recover(data_root+"model/",data_root+"model/mesh.elastic","./tempt/test5",
		  W,lambda,1e-3,1e-3,1e-3,1e16);

}

int main(int argc, char *argv[]){

  // recoverCoarse();
  recoverFine();
}
