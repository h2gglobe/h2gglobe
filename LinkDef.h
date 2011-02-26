#ifndef ROOT_TREE_VECTOR_LINKDEF_H 
#define ROOT_TREE_VECTOR_LINKDEF_H 1

#ifdef __CINT__

//#pragma link off all classes;

#pragma link C++ class std::vector<std::vector<unsigned short> >+;
#pragma link C++ class std::vector<std::vector<int> >+;
#pragma link C++ class std::vector<unsigned short>+;
#pragma link C++ class std::vector<std::string>+;
//#pragma link C++ class std::map<std::string, int>+;

#endif

#endif // ROOT_TREE_VECTOR_LINKDEF_H
