#include <ElasticForceTetFullStVK.h>
using namespace UTILITY;

void ElasticForceTetFullStVK::force_tet(mat3x4& f,const int& i,const VectorXd& X){

  Matrix3d  F;
  _def_grad.evalFe(F,X,i);
  Matrix3d  P;
  PF(P,F,_vol_mesh->material()._G[i],_vol_mesh->material()._lambda[i]);
  f.block<3,3>(0,0)=-_volume[i]*P*_def_grad.invDm()[i].transpose();
  f.col(3)=-(f.col(0)+f.col(1)+f.col(2));
}

void ElasticForceTetFullStVK::forceDerivX_tet(TetDF &df, const int& i, const VectorXd& X){

	Matrix3d  F;
	_def_grad.evalFe(F,X,i);

	const derivF& dF=_def_grad.dF()[i];
	for(int x=0;x<4;x++)
	for(int d=0;d<3;d++){

	  const Matrix3d& dFxd=dF._dF[x][d];
	  Matrix3d  deriv;
	  dPF(deriv,dFxd,F,_vol_mesh->material()._G[i],_vol_mesh->material()._lambda[i]);
	  deriv=-_volume[i]*deriv*_def_grad.invDm()[i].transpose();

	  for(int fx=0;fx<3;fx++)
		for(int fd=0;fd<3;fd++)
		  df.df[fx][x](fd,d)=deriv(fd,fx);
	}

	df.df[3][0]=-df.df[0][0]-df.df[1][0]-df.df[2][0];
	df.df[3][1]=-df.df[0][1]-df.df[1][1]-df.df[2][1];
	df.df[3][2]=-df.df[0][2]-df.df[1][2]-df.df[2][2];
	df.df[3][3]=-df.df[0][3]-df.df[1][3]-df.df[2][3];
}

void ElasticForceTetFullStVK::forceDerivXdX_tet(mat3x4& dfdX,const int& i,const VectorXd& dx,const VectorXd& X){

  Matrix3d  F;
  _def_grad.evalFe(F,X,i);

  Matrix3d  dF;
  _def_grad.evalFe(dF,dx,i);
		
  Matrix3d  deriv;
  dPF(deriv,dF,F,_vol_mesh->material()._G[i],_vol_mesh->material()._lambda[i]);
		
  dfdX.block<3,3>(0,0)=-_volume[i]*deriv*_def_grad.invDm()[i];
  dfdX.col(3)=-(dfdX.col(0)+dfdX.col(1)+dfdX.col(2));
}

void ElasticForceTetFullStVK::W(double& W,const Matrix3d& F,const double& G,const double& lambda,const double& V) const{

  const Matrix3d  E=0.5f*(F.transpose()*F-Matrix3d::Identity());
  const double tr=E.trace();
  W=G*E.squaredNorm()+lambda*0.5f*tr*tr;
  W*=V;
}

void ElasticForceTetFullStVK::PF(Matrix3d& P,const Matrix3d& F,const double& G,const double& lambda) const{

  const Matrix3d  E=(F.transpose()*F-Matrix3d::Identity())*0.5f;
  P=F*E*(2.0f*G)+(E.trace()*lambda)*F;
}

void ElasticForceTetFullStVK::dPF(Matrix3d& deriv,const Matrix3d& dF,const Matrix3d& F,const double& G,const double& lambda) const{

  const Matrix3d  E=(F.transpose()*F-Matrix3d::Identity())*0.5f;
  const Matrix3d  derivE=(dF.transpose()*F+F.transpose()*dF)*0.5f;
  deriv=2.0f*G*(dF*E+F*derivE)+lambda*(E.trace()*dF+derivE.trace()*F);
}
