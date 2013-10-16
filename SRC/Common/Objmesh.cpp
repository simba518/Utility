#include <boost/filesystem.hpp>
#include <Objmesh.h>
#include <AuxTools.h>
using namespace boost::filesystem;
using namespace UTILITY;

void ObjMtl::setDefault(){

  name = "defualt";

  ambient[0] = 0.0f;       ambient[0] = 0.1f;       ambient[0] = 0.0f;
  diffuse[0] = 0.0f;       diffuse[0] = 0.6f;       diffuse[0] = 0.0f;
  specular[0] = 0.35f;      specular[0] = 0.35f;      specular[0] = 0.35f;
  transmittance[0] = 1.0f; transmittance[0] = 1.0f; transmittance[0] = 1.0f;
  emission[0] = 0.0f;      emission[0] = 0.0f;      emission[0] = 0.0f;

  shininess = 200;
  ior = 1;
}

// v : 0
// vn: 1
// vt: 2
// f : 3
// other: -1
int Objmesh::type(const string &line)const{
  
  int t = -1;
  if(line.size() > 2){
	if('v' == line[0]){
	  t = 0;
	  if('n' == line[1])
		t = 1;
	  else if('t' == line[1])
		t = 2;
	}else if('f' == line[0])
	  t = 3;
  }
  return t;
}

bool Objmesh::getFaceVert(string &line,int &v1,int &v2,int &v3,int &vn1,int &vn2,int &vn3)const{

  vn1=-1, vn2=-1, vn3=-1;
  int numSlash = 0;
  for (int i = 0; i < line.size(); ++i){
    if('/' == line[i]){
	  line[i] = ' ';
	  numSlash ++;
	  if('/' == line[i+1])
		line[i+1] = ' ';
	}
  }
  istringstream ss(line);
  string temp;
  ss >> temp;
  if (numSlash == 0){
	ss >> v1 >> v2 >> v3;
  }else if (numSlash == 3){
	ss >> v1>> vn1 >> v2>> vn2 >> v3 >> vn3;
  }else if (numSlash == 6){
	ss >> v1 >> temp >> vn1 >> v2 >> temp >> vn2  >> v3 >> temp >> vn3;
  }else{
	ERROR_LOG("the format of the face is not supported, only triangle mesh supported.");
	return false;
  }
  v1 -= 1;
  v2 -= 1;
  v3 -= 1;
  vn1 -= 1;
  vn2 -= 1;
  vn3 -= 1;

  return true;
}

// load triangle obj file
// 1. support only triangle mesh.
// 2. load only one mtl.
// 3. no textures.
bool Objmesh::load(const string fname){

  std::ifstream inf;
  inf.open(string(fname).c_str());
  if(!inf.is_open()){
	ERROR_LOG("failed to open file for output: "<<fname); 
	return false;  
  }
  
  // count the number of v,vn,vt and f.
  int numType[4] = {0}; // num of v,vn,vt,f
  string line;
  string mtlfile;
  while(!inf.eof()&&getline(inf,line)){
	const int t = type(line);
	if(t >=0 && t < 4){
	  numType[t]++;
	}else if(line.find("mtllib") != std::string::npos){
	  const string dir = string(path(fname).parent_path().c_str())+"/";
	  mtlfile = dir + line.substr(7,line.size()-7);
	}
  }
  
  _verts.resize(numType[0]*3);  
  _verts.setZero();
  _vertNormal.resize(numType[1]*3);
  _vertNormal.setZero();
  _faces.resize(numType[3]*3);
  _faces.setZero();

  // load real data
  inf.clear();
  inf.seekg(0,ios::beg);
  int v=0,vn=0,f=0;
  vector<int> normalIndex;
  normalIndex.reserve(numType[3]*3);
  bool succ = true;
  while(!inf.eof()&&getline(inf,line) && succ){
	const int t = type(line);
	if(t >= 0){
	  istringstream ss(line);
	  string temp;
	  ss >> temp;
	  switch (t){
	  case 0:{
		ss >> _verts[v] >> _verts[v+1] >> _verts[v+2];
		v += 3;
		break;
	  };
	  case 1:{
		ss >> _vertNormal[vn] >> _vertNormal[vn+1] >> _vertNormal[vn+2];
		vn += 3;
		break;
	  };
	  case 3:{
		int vn1,vn2,vn3;
		succ = getFaceVert(line,_faces[f],_faces[f+1],_faces[f+2],vn1,vn2,vn3);
		f += 3;
		if(succ && vn1 >= 0){
		  normalIndex.push_back(vn1);
		  normalIndex.push_back(vn2);
		  normalIndex.push_back(vn3);
		}
		break;
	  };
	  }
	}
  }

  // average the normals.
  if(normalIndex.size() == _faces.size()){

	vector<float> counter(_verts.size()/3,0);
	const Eigen::VectorXd normal = _vertNormal;
	_vertNormal.resize(_verts.size());
	_vertNormal.setZero();
	for (int i = 0; i < normalIndex.size(); ++i){
	  const int vn = normalIndex[i];
	  const int vi = _faces[i];
	  assert_in(vn*3,0,normal.size()-3);
	  assert_in(vi*3,0,_verts.size()-3);
	  _vertNormal.segment(vi*3,3) += normal.segment(vn*3,3);
	  counter[vi] ++;
	}
	for (int i = 0; i < counter.size(); ++i){
	  assert_gt(counter[i],0.0f);
	  _vertNormal.segment(i*3,3) /= counter[i];
	}
  }

  inf.close();

  // load mtl
  if(mtlfile.size() > 0){
	if(!loadMtl(mtlfile)){
	  _mtl.setDefault();
	  WARN_LOG("failed to load the mtl from "<<mtlfile<<"default mtl will be used.");
	}
  }

  return succ;
}

bool Objmesh::write(const string fname)const{

  ofstream outf;
  outf.open(fname.c_str());
  if(!outf.is_open()){
	ERROR_LOG("failed to open file "<< fname <<" for writing.");
	return false;  
  }

  const int nv = getVertsNum();
  const int nf = getFacesNum();
  outf << "# number of vertices "<< nv <<endl;
  outf << "# number of faces "<< nf<< endl;
  if(_vertNormal.size() == nv*3)
	outf << "# number of normals " << nv << endl;
  outf << endl;

  for (int i = 0; i < nv; ++i){
	outf<< "v " << _verts[i*3] << " ";
	outf<< _verts[i*3+1] << " ";
	outf<< _verts[i*3+2] << endl;
  }
  outf << endl;

  if(_vertNormal.size() == nv*3){
	for (int i = 0; i < nv; ++i){
	  outf<< "vn " << _vertNormal[i*3] << " ";
	  outf<< _vertNormal[i*3+1] << " ";
	  outf<< _vertNormal[i*3+2] << endl;
	}
	outf << endl;
	for (int i = 0; i < nf; ++i){
	  outf<<"f "<< _faces[i*3]+1<<"//"<<_faces[i*3]+1<<" ";
	  outf<< _faces[i*3+1]+1<<"//"<<_faces[i*3+1]+1<<" ";
	  outf<< _faces[i*3+2]+1<<"//"<<_faces[i*3+2]+1<<endl;
	}
  }else{
	for (int i = 0; i < nf; ++i){
	  outf<<"f "<< _faces[i*3]+1<<" ";
	  outf<< _faces[i*3+1]+1<<" ";
	  outf<< _faces[i*3+2]+1<<endl;
	}
  }

  outf.close();
  return true;
}

bool Objmesh::loadMtl(const string fname){
  
  std::ifstream inf;
  inf.open(string(fname).c_str());
  if(!inf.is_open()){
	ERROR_LOG("failed to open file for output: "<<fname);
	return false;  
  }
  string key;
  while(!inf.eof()&&(inf >> key).good()){
	if(key.size()<=0||'#' == key[0]){
	  continue;
	}else if(string("Kd") == key){
	  inf >> _mtl.diffuse[0] >> _mtl.diffuse[1] >> _mtl.diffuse[2];
	}else if(string("Ka") == key){
	  inf >> _mtl.ambient[0] >> _mtl.ambient[1] >> _mtl.ambient[2];
	}else if(string("Tf") == key){
	  inf >> _mtl.transmittance[0] >> _mtl.transmittance[1] >> _mtl.transmittance[2];
	}else if(string("Ks") == key){
	  inf >> _mtl.specular[0] >> _mtl.specular[1] >> _mtl.specular[2];
	}else if(string("Ke") == key){
	  inf >> _mtl.emission[0] >> _mtl.emission[1] >> _mtl.emission[2];
	}else if(string("Ns") == key){
	  inf >> _mtl.shininess;
	}else if(string("Ni") == key){
	  inf >> _mtl.ior;
	}else if(string("newmtl") == key){
	  inf >> _mtl.name;
	}
  }
  inf.close();
  return true;
}
