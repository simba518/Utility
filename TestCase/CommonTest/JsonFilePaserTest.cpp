#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <JsonFilePaser.h>

struct JsonFilePaserTestInit{
  JsonFilePaserTestInit(){
	jsonFilePath = "./TestCase/TestData/json.in";
	ASSERT_MSG("failed to open file:"<<jsonFilePath,jsonFile.open(jsonFilePath));
  }
  std::string jsonFilePath;
  UTILITY::JsonFilePaser jsonFile;
};

BOOST_AUTO_TEST_SUITE(JsonFilePaser)

BOOST_FIXTURE_TEST_CASE(readFloat, JsonFilePaserTestInit)
{
  float f;
  ASSERT(jsonFile.read("floatnumber",f));
  ASSERT_EQ(f,100.0);
}

BOOST_FIXTURE_TEST_CASE(readVector, JsonFilePaserTestInit)
{
  std::vector<int> v, cv;
  cv.push_back(1);
  cv.push_back(2);
  cv.push_back(3);
  ASSERT(jsonFile.read("intArray",v));
  ASSERT_EQ(v.size(),cv.size());
  ASSERT_EQ_SMALL_VEC(v,cv,v.size());
}

BOOST_FIXTURE_TEST_CASE(readFilePath, JsonFilePaserTestInit)
{
  std::string filename;
  ASSERT(jsonFile.readFilePath("filename",filename,false));
  ASSERT_NE(filename,std::string("f.txt"));
  ASSERT_EQ(filename,jsonFile.getFileDir()+std::string("f.txt"));
}

BOOST_AUTO_TEST_SUITE_END()
