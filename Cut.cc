#include "Cut.h"
#include <iostream>

Cut::Cut() {}

Cut::~Cut() {}

void Cut::Print() {
  std::cout<<"cutStruct test print"<<std::endl;
  std::cout<<"name "<<name<<std::endl;
  std::cout<<"fromright "<<fromright<<std::endl;
  std::cout<<"finalcut "<<finalcut<<std::endl;
  std::cout<<"index "<<index<<std::endl;
  std::cout<<"ncat "<<ncat<<std::endl;
  for (int j=0; j<ncat; j++) {
    if(fromright==2) {
      std::cout<<j<<" cut range "<<cutintervall[j]<<" "<<cutintervalh[j]<<std::endl;
    } else {
      std::cout<<j<<" cut "<<cut[j]<<std::endl;
    }
  }
}
