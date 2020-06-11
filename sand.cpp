// basic file operations
#include <iostream>
#include <fstream>
using namespace std;

void printFunction(int k){

  double result = (double) 1/k;
  cout << result << endl;
  
}

int main () {
  
  int k = 2;
  printFunction(k);
  return 0;
  
}
