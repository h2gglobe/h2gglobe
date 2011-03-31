#include "CounterContainer.h"
#include <utility>
#include <iostream>

CounterContainer::CounterContainer(int n) {
  histVal = n;
}

CounterContainer::~CounterContainer() 
{}

void CounterContainer::Add(std::string name, int categories, 
			   std::string denom1, std::string denom2, std::string denom3) {
  int d1=-1, d2=-1, d3=-1;

  std::map<std::string, std::vector<int> >::iterator it;
  it = c.find(denom1);
  if (it != c.end())
    d1 = std::distance(c.begin(), it);

  it = c.find(denom2);
  if (it != c.end())
    d2 = std::distance(c.begin(), it);

  it = c.find(denom3);
  if (it != c.end())
    d3 = std::distance(c.begin(), it);

  std::vector<int> temp;
  for (unsigned int i=0; i<categories; i++)
    temp.push_back(0);

  c.insert(std::pair<std::string, std::vector<int> >(name, temp));
  denom1_.push_back(d1);
  denom2_.push_back(d2);
  denom3_.push_back(d3);
}

void CounterContainer::Fill(std::string name, int category) {
  Fill(name, category,1.0);
}

void CounterContainer::Fill(std::string name, int category, float weight) {

  std::map<std::string, std::vector<int> >::iterator it = c.find(name);
  
  if (it != c.end()) {
    c[name][category]++;
  }
}

unsigned int CounterContainer::ncat(int length) {
  std::map<std::string, std::vector<int> >::const_iterator it = c.begin();
  
  for (int i=0; i<length; i++)
    it++;

  return it->second.size();
}

std::vector<int> CounterContainer::operator[](unsigned int length) {

  std::map<std::string, std::vector<int> >::const_iterator it = c.begin();
  
  for (int i=0; i<length; i++)
    it++;

  return it->second;
}

std::string CounterContainer::name(int length) {

  std::map<std::string, std::vector<int> >::const_iterator it = c.begin();
  for (int i=0; i<length; i++)
    it++;

  return it->first;
}
