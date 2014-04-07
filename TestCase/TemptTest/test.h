
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


#include <iostream>
#include <cmath>
#include <omp.h>
#include <stdlib.h>
#include <float.h>
#include "QPSolver.h"
using namespace std;

float ScalarUtil<float>::scalar_max=FLT_MAX;
float ScalarUtil<float>::scalar_eps=1E-5f;

double ScalarUtil<double>::scalar_max=DBL_MAX;
double ScalarUtil<double>::scalar_eps=1E-9;

int main(int argc, char *argv[]){

  FixedSparseMatrix<double> A;
  MPRGPQPSolver<double>::Vec B(1),L(1),H(1);
  A.resize(1,1);
  A.insert(0,0) = 1;
  B << -1;
  L << 1.0f;
  H << 10;

  // solve for min 1/2*x^t*A*x-x^t*B, s.t.  L<=x<=H.
  MPRGPQPSolver<double> f(A,B,L,H);
  f.setCallback(boost::shared_ptr<Callback<double> >(new Callback<double>()));
  MPRGPQPSolver<double>::Vec x(1);
  x << 3.0;
  f.solve(x);
  printf("\n");
  cout << x.transpose() << endl;

  return 0;
}
