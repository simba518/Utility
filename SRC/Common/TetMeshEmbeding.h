#ifndef _TETMESHEMBEDING_H_
#define _TETMESHEMBEDING_H_

#include <TetMesh.h>
#include <ObjMesh.h>

namespace UTILITY{

  class TetMeshEmbeding{

  public:
	TetMeshEmbeding(){}
	TetMeshEmbeding(pTetMesh_const tetMesh,pObjMesh objMesh){
	  setTetMesh(tetMesh);
	  setObjMesh(objMesh);
	}
	void setTetMesh(pTetMesh_const tetMesh){
	  assert(_tetMesh);
	  _tetMesh = tetMesh;
	}
	void setObjMesh(pObjMesh objMesh){
	  assert(_objMesh);
	  _objMesh = objMesh;
	  _objRestVerts = _objMesh->getVerts();
	}

	void buildInterpWeights(){
	  assert(_tetMesh);
	  assert(_objMesh);
	  _tetMesh->buildInterpWeights(_objMesh->getVerts(),_nodes,_weights);
	}
	void interpolate(const VectorXd& u){
	  assert(_tetMesh);
	  assert(_objMesh);
	  _uTarget.resize(_objRestVerts.size());
	  _tetMesh->interpolate(_nodes,_weights,u,_uTarget);
	  _uTarget += _objRestVerts;
	  _objMesh->setVerts(_uTarget);
	}

  private:
	pTetMesh_const _tetMesh;
	pObjMesh _objMesh;
	VectorXd _objRestVerts;

	vector<int> _nodes;
	VectorXd _weights;

	VectorXd _uTarget;
  };
  
  typedef boost::shared_ptr<TetMeshEmbeding> pTetMeshEmbeding;

}//end of namespace

#endif /* _TETMESHEMBEDING_H_ */
