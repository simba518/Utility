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

void recover(MaterialFitting *mtlfit_k,
			 const string data_root, const string mtl_file, const string save_to,
			 const MatrixXd &W, const VectorXd &lambda, 
			 const double Gs=1e-3, const double Gl=1e-3,
			 const double lower=1, const double upper=1e10, const bool use_ipopt=true){

  { // fit G, l
  	INFO_LOG("fit G, l");
  	mtlfit_k->loadTetMesh(data_root+"mesh.abq");
  	mtlfit_k->loadMtl(mtl_file);
  	mtlfit_k->loadFixednodes(data_root+"/con_nodes.bou");
  	mtlfit_k->setWLambda(W, lambda);
  	mtlfit_k->computeK();
  	mtlfit_k->computeM();
  	mtlfit_k->removeFixedDOFs();

  	mtlfit_k->setBounds(lower,upper);
  	mtlfit_k->setMuSmoothGL(Gs, Gl);
  	mtlfit_k->setMuSmoothEv(0, 0);
  	mtlfit_k->setMuSmoothDensity(0.0f);
  	mtlfit_k->setMuAverageDensity(0.0f);
	// mtlfit_k->scale();
  	mtlfit_k->assembleObjfun();

	if (use_ipopt){
	  mtlfit_k->solveByIpopt();
	  // mtlfit_k->unScale();
	  mtlfit_k->saveResults(save_to);
	}else{
	  mtlfit_k->solveByMPRGP(save_to+"_x.b");
	  // mtlfit_k->unScale();
	  mtlfit_k->saveResults(save_to+"_MPRGP_");
	}
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

void loadAndScale(const string data_root_opt, MatrixXd &W, VectorXd &lambda,const double h, const double rho){
  
  MatrixXd S, eig_W;
  VectorXd eig_lambda;
  bool succ = load(data_root_opt+"opt_scaled_W.b",eig_W); assert(succ);
  succ = load(data_root_opt+"opt_S.b",S); assert(succ);
  const MatrixXd t = eig_W;
  eig_W = t.leftCols(S.rows())*S;
  eig_W *= 1.0f/rho;

  for (int i = 0; i < eig_W.cols(); ++i){
    cout << eig_W.col(i).norm() << endl;
  }

  succ = load(data_root_opt+"opt_lambda.b",eig_lambda); assert(succ);
  eig_lambda *= rho/(h*h);

  MatrixXd Z,F;
  succ = load(data_root_opt+"opt_Z.b",Z); assert(succ);
  succ = load(data_root_opt+"opt_F.b",F); assert(succ);

  const int T = 150;
  double s_max = -1.0f;
  for (int i = 0; i < Z.rows(); ++i){
	const double dz = Z.row(i).head(T).norm();
	const double df = F.row(i).head(T).norm();
	double s = dz/(df*eig_lambda[i]);
	s = s*s;
	eig_W.col(i) *= s;
	if (s > s_max)
	  s_max = s;
	cout<< i << ",\t" << dz<<",\t" << df << ",\t" << s << ",\t"<< eig_lambda[i] << "\n";
  }
  assert_gt(s_max, 0.0f);
  cout << "smax: " << s_max << "\n\n";
  eig_W *= 1.0f/s_max;

  W = eig_W;
  lambda = eig_lambda;
  cout<< "\n\nLambda: " << lambda.transpose() << "\n\n";
}

void recoverCoarseNonSmooth(){

  TRACE_FUN();

  const string data_root = "/home/simba/Workspace/SolidSimulator/data/beam-coarse/model/";
  const string data_root_opt ="/home/simba/Workspace/AnimationEditor/Data/beam-coarse/model/";

  MatrixXd W, eig_W;
  VectorXd lambda, eig_lambda;
  loadAndScale(data_root_opt, eig_W, eig_lambda, 0.03f,1.0f);

  { //save real mateial
	TetMesh tetmesh;
	bool succ = tetmesh.load(data_root+"mesh.abq"); assert(succ);
	succ = tetmesh.loadElasticMtl(data_root+"/mesh.elastic"); assert(succ);
	succ = tetmesh.writeElasticMtlVTK("./tempt/coarse/coarse_real"); assert(succ);
  }

  { // load real value
  	// bool succ = load(data_root+"/tempt_eigenvectors.b",eig_W); assert(succ);
  	// succ = load(data_root+"/tempt_eigenvalues.b",eig_lambda); assert(succ);
  	// cout << "lambda: " << eig_lambda.transpose() << "\n\n";
  }
  
  W = eig_W;
  lambda = eig_lambda;
  VectorXd mu_e(3);
  mu_e << 1e-2, 1e-6, 1e-4;
  for (int i = 0; i < mu_e.size(); ++i){
	MaterialFitting_MA_K mtlfit_k;
	mtlfit_k.useHessian(true);
	recover(&mtlfit_k,data_root, data_root_opt+"/mesh.elastic", 
			"./tempt/coarse/test_coarse_non_smooth"+TOSTR(mu_e[i]), 
			W, lambda,mu_e[i],mu_e[i],0,1e16,true);
  }
}

void recoverFineSmooth(){

  TRACE_FUN();

  const string data_root = "/home/simba/Workspace/SolidSimulator/data/beam/model/";
  const string data_root_opt ="/home/simba/Workspace/AnimationEditor/Data/beam_fine/model/";

  MatrixXd eig_W, W;
  VectorXd eig_lambda, lambda;
  loadAndScale(data_root_opt+"smooth/", eig_W, eig_lambda, 0.04f, 100.0f);

  W = eig_W;
  lambda = eig_lambda;
  VectorXd mu_e(5);
  mu_e << 1e-5,1e-6,1e-8,1e-3,1e-4;
  for (int i = 0; i < mu_e.size(); ++i){

	MaterialFitting_MA_K mtlfit_k;
	mtlfit_k.useHessian(false);
	recover(&mtlfit_k,data_root, data_root+"/tempt_mesh.elastic", 
			"./tempt/fine_smooth/test_fine_smooth"+TOSTR(mu_e[i]), 
			W, lambda,mu_e[i],mu_e[i],10.0,1e12,true);
  }
}

void recoverFineNonSmooth(){

  TRACE_FUN();

  const string data_root = "/home/simba/Workspace/SolidSimulator/data/beam/";
  const string data_root_opt ="/home/simba/Workspace/AnimationEditor/Data/beam_fine/model/";

  MatrixXd W, eig_W;
  VectorXd lambda, eig_lambda;
  loadAndScale(data_root_opt, eig_W, eig_lambda, 0.04f,100.0f);

  { //save real mateial
	TetMesh tetmesh;
	bool succ = tetmesh.load(data_root+"/model/mesh.abq"); assert(succ);
	succ = tetmesh.loadElasticMtl(data_root+"/model/mesh.elastic"); assert(succ);
	succ = tetmesh.writeElasticMtlVTK("./tempt/fine/coarse_real"); assert(succ);
  }

  { // load real value
  	bool succ = load(data_root+"tempt/eigenvectors.b",eig_W); assert(succ);
  	succ = load(data_root+"tempt/eigenvalues.b",eig_lambda); assert(succ);
  	cout << "lambda: " << lambda.transpose() << "\n\n";
  }

  W = eig_W.leftCols(10);
  lambda = eig_lambda.head(10);
  VectorXd mu_e(3);
  mu_e << 1e-3, 1e-2, 1e-4;
  for (int i = 0; i < mu_e.size(); ++i){
	
	MaterialFitting_MA_K mtlfit_k;
	mtlfit_k.useHessian(false);
	recover(&mtlfit_k,data_root+"/model/", data_root+"model/tempt_mesh.elastic", 
			"./tempt/fine/test_fine_non_smooth"+TOSTR(mu_e[i]), 
			W, lambda,mu_e[i],mu_e[i],100,1e12,true);
  }

}

int main(int argc, char *argv[]){

  // recoverCoarseNonSmooth();
  recoverFineNonSmooth();
  // recoverFineSmooth();
}
