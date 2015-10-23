#include <iostream>
#include "../../treader.h"

using namespace std;
using namespace omd;

int main() {
  TableReader tab("gcu.dat");
  cout<<"# range of table ("<<tab.min_range()<<","<<tab.max_range()<<")\n";
  
  double step=(tab.max_range()-tab.min_range())/1000;
  
  for(double x=tab.min_range();x<tab.max_range();x+=step) {
    cout<<x<<" "<<tab.read(x)<<" "<<tab.dread(x)<<" "<<tab.dread2(x)<<"\n";
  }
  
  return 0;
}
