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

void MaterialFitting::assembleObjfun(){
  
  TRACE_FUN();
  objfun = 0.0f;
  addSmoothObjfun(objfun);

  const SXMatrix M1 = K.mul(W)-M.mul(W.mul(lambda));
  for (int i = 0; i < M1.size1(); ++i){
    for (int j = 0; j < M1.size2(); ++j)
	  objfun += M1.elem(i,j)*M1.elem(i,j);
  }

  objfun *= 0.5f;
}

void MaterialFitting::solve(){
  
  TRACE_FUN();
  INFO_LOG("init solver");
  VSX G_Lame_rho;
  initAllVariables(G_Lame_rho);
  CasADi::SXFunction fun = CasADi::SXFunction(G_Lame_rho,objfun);
  CasADi::IpoptSolver solver = CasADi::IpoptSolver(fun);

  solver.setOption("generate_hessian",use_hessian);
  solver.setOption("tol",1e-18);
  solver.setOption("max_iter",1000);

  Timer timer;
  timer.start();
  solver.init();
  timer.stop("solver.init() ");

  VectorXd init_x;
  init_x.resize(G_Lame_rho.size());
  const vector<double> &g = tetmesh->material()._G;
  const vector<double> &la = tetmesh->material()._lambda;
  const vector<double> &rh = tetmesh->material()._rho;
  const int num_tet = tetmesh->tets().size();
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
  for (int i = 0; i < lower.size(); ++i){
    lower[i] = 0.0f;
  }
  solver.setInput(lower,CasADi::NLP_LBX);

  // solving
  INFO_LOG("solving");
  timer.start();
  solver.solve();
  timer.stop("solver.solve() ");
  rlst.resize(G_Lame_rho.size());
  solver.getOutput(rlst,CasADi::NLP_X_OPT);
}

void MaterialFitting::saveResults(const string filename)const{
  
  // convert to Young's(E) and Poisson(v)
  const ElasticMaterial<double> mtl_0 = tetmesh->material();
  const int num_tet = tetmesh->tets().size();
  const vector<double> vx = rlst;

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

  bool succ = tetmesh->writeElasticMtlVTK(filename);
  assert(succ);

  tetmesh->material() = mtl_0;
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
		objfun += 0.25f*mu_neigh*(vi+vj)*(G[i]-G[j])*(G[i]-G[j]);
		objfun += 0.25f*mu_neigh*(vi+vj)*(Lame[i]-Lame[j])*(Lame[i]-Lame[j]);
		objfun += 0.25f*mu_neigh*(vi+vj)*(rho[i]-rho[j])*(rho[i]-rho[j]);
	  }
	}
  }
}

void MaterialFittingM2::assembleObjfun(){
  
  TRACE_FUN();
  objfun = 0.0f;
  addSmoothObjfun(objfun);

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
