#ifndef _TETMESHEMBEDING_H_
#define _TETMESHEMBEDING_H_

#include <stdio.h>
#include <TetMesh.h>
#include <Objmesh.h>
#include <ObjFileIO.h>

namespace UTILITY{

  class TetMeshEmbeding{

  public:
	TetMeshEmbeding(){
	  _objMesh = pObjmesh(new Objmesh());
	  _tetMesh = pTetMesh(new TetMesh());
	}
	TetMeshEmbeding(pTetMesh tetMesh,pObjmesh objMesh){
	  setTetMesh(tetMesh);
	  setObjmesh(objMesh);
	}
	void setTetMesh(pTetMesh tetMesh){
	  assert(_tetMesh);
	  _tetMesh = tetMesh;
	}
	void setObjmesh(pObjmesh objMesh){
	  assert(_objMesh);
	  _objMesh = objMesh;
	  _objRestVerts = _objMesh->getVerts();
	}

	pTetMesh getTetMesh(){
	  return _tetMesh;
	}
	pObjmesh getObjMesh(){
	  return _objMesh;
	}
	pTetMesh_const getTetMesh()const{
	  return _tetMesh;
	}
	pObjmesh_const getObjMesh()const{
	  return _objMesh;
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

	bool loadObjMesh(const string filename){
	  assert(_objMesh);
	  return load(filename,*_objMesh);
	}
	bool loadTetMesh(const string filename){
	  assert(_tetMesh);
	  return _tetMesh->load(filename);
	}
	
	bool loadWeights(const string fname){

	  FILE *fin = fopen(fname.c_str(),"ra");
	  if (!fin){
		ERROR_LOG("unable to open file"<<fname);
		return false;
	  }
	  
	  const int numElementVertices = 4;
	  const int numTargetLocations = _objMesh->getVertsNum();
	  _nodes.resize(numElementVertices*numTargetLocations);
	  _weights.resize(_nodes.size());

	  int numVertices = -1;
	  int currentVertex;

	  // read the elements one by one and accumulate entries
	  while (numVertices < numTargetLocations-1){
		numVertices++;
		if ( feof(fin) ){
		  ERROR_LOG("interpolation file is too short.");
		  return false;
		}
		fscanf(fin, "%d", &currentVertex);
		if (currentVertex != numVertices){
		  ERROR_LOG("consecutive vertex index at position: "<<currentVertex<<" mismatch.");
		  return false;
		}
		for(int j=0; j<numElementVertices; j++)
		  fscanf(fin,"%d %lf", &(_nodes[currentVertex * numElementVertices+j]),
				 &(_weights[currentVertex*numElementVertices+j]));
		fscanf(fin,"\n");
	  }

	  fclose(fin);
	  return true;
	}
	bool writeWeights(const string fname)const{
	  return false;
	}

	void getBBox(double min[3],double max[3])const{
	  assert(_tetMesh);
	  assert(_objMesh);
	  BBoxD b1 = _tetMesh->getBBox();
	  const BBoxD b2 = _objMesh->getBBox();
	  b1.add(b2);
	  b1.getMaxConner(max);
	  b1.getMinConner(min);
	}
	double getMaxRadius()const{
	  double min[3],max[3];
	  double r = max[0] - min[0];
	  r = (r>=(max[1]-min[1]) ? r:(max[1]-min[1]));
	  r = (r>=(max[2]-min[2]) ? r:(max[2]-min[2]));
	  return r >= 0 ? r : 0;
	}

	const vector<int> &getInterpNodes()const{
	  return _nodes;
	}
	const VectorXd &getInterpWeights()const{
	  return _weights;
	}

  private:
	pTetMesh _tetMesh;
	pObjmesh _objMesh;
	VectorXd _objRestVerts;

	vector<int> _nodes;
	VectorXd _weights;

	VectorXd _uTarget;
  };
  
  typedef boost::shared_ptr<TetMeshEmbeding> pTetMeshEmbeding;
  typedef boost::shared_ptr<const TetMeshEmbeding> pTetMeshEmbeding_const;

}//end of namespace

#endif /* _TETMESHEMBEDING_H_ */
