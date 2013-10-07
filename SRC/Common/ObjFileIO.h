//
// Copyright 2012-2013, Syoyo Fujita.
//
// Licensed under 2-clause BSD liecense.
//

#ifndef _OBJFILEIO_H_
#define _OBJFILEIO_H_

#include <string>
#include <vector>
#include <map>
#include <Objmesh.h>
#include <Log.h>

namespace UTILITY {

  typedef ObjMtl material_t;

  typedef struct {
    std::vector<float>          positions;
    std::vector<float>          normals;
    std::vector<float>          texcoords;
    std::vector<unsigned int>   indices;
  } mesh_t;

  typedef struct {
    std::string  name;
    material_t   material;
    mesh_t       mesh;
  } shape_t;

  /// Loads .obj from a file.
  /// 'shapes' will be filled with parsed shape data
  /// The function returns error string.
  /// Returns empty string when loading .obj success.
  /// 'mtl_basepath' is optional, and used for base path for .mtl file.
  std::string LoadObj( std::vector<shape_t>& shapes,   // [output]
					   const char* filename,
					   const char* mtl_basepath = NULL);

  // mesh io
  bool load(const string fname,Objmesh &mesh);

  inline bool write(const string fname,const Objmesh &mesh){
	ERROR_LOG("undefined function.");
	return false;
  }

};

#endif /* _OBJFILEIO_H_ */
