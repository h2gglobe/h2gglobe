#include "CounterContainer.h"
#include <utility>
#include <iostream>
#include <assert.h>

CounterContainer::CounterContainer(int n) {
  histVal = n;
}

CounterContainer::~CounterContainer() 
{}

void CounterContainer::Add(std::string name, int categories, 
			   std::string denom1, std::string denom2, std::string denom3) {
  std::string d1(""), d2(""), d3("");
  
  std::vector<float> temp;
  for (unsigned int i=0; i<categories; i++)
    temp.push_back(0.0);

  c.push_back(temp);
  names.push_back(name);
  
  std::vector<std::string> stemp;
  stemp.push_back(d1);
  stemp.push_back(d2);
  stemp.push_back(d3);
  denoms_.push_back(stemp);
}

void CounterContainer::Fill(std::string name, int category) {
  Fill(name, category, 1.0);
}

void CounterContainer::Fill(std::string name, int category, float weight) {
  
  int index = -1;
  for (unsigned int i=0; i<names.size(); i++) {
    if (names[i] == name) {
      index = i;
      break;
    } 
  }
  
  if (index != -1)
    c[index][category] += weight;
  else
    std::cout << "Wrong counter name: " <<name<< std::endl;
}

unsigned int CounterContainer::ncat(unsigned int length) {
  return c[length].size();
}

float CounterContainer::tot(unsigned int length) {
  std::vector<float> temp = c[length];

  float t = 0;
  for(unsigned int i=0; i<temp.size(); i++)
    t += temp[i];

  return t;
}

std::vector<float> CounterContainer::operator[](unsigned int length) {

  return c[length];
}

std::string CounterContainer::denomName(unsigned int length, unsigned int den) {

  return denoms_[length][den];
}

std::string CounterContainer::name(unsigned int length) {

  return names[length];
}

float CounterContainer::efficiency(unsigned int index, unsigned int cat, unsigned int denom_type) {

  if (index == -1)
    return -1.;  
  
  if (denom_type > 2 || denom_type < 0) {
    std::cout << "Wrong denominator" << std::endl;
    return -1;
  }
  
  if (denoms_[index][denom_type] == "")
    return -1.;
  
  int den_index = -1;
  for (unsigned int i=0; i<names.size(); i++) {
    if (names[i] == denoms_[index][denom_type]) {
      den_index = i;
      break;
    }
  }

  if (den_index == -1) {
    std::cout << "Wrong denominator" << std::endl;
    return -1;
  } else
    return c[index][cat]/c[den_index][cat];
}


float CounterContainer::efficiency(unsigned int index, unsigned int denom_type) {

  if (index == -1)
    return -1.;  
  
  if (denom_type > 2 || denom_type < 0) {
    std::cout << "Wrong denominator" << std::endl;
    return -1;
  }
  
  if (denoms_[index][denom_type] == "")
    return -1.;
  
  int den_index = -1;
  for (unsigned int i=0; i<names.size(); i++) {
    if (names[i] == denoms_[index][denom_type]) {
      den_index = i;
      break;
    }
  }

  if (den_index == -1) {
    std::cout << "Wrong denominator" << std::endl;
    return -1;
  } else
    return tot(index)/tot(den_index);
}
