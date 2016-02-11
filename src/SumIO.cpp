// IO functions of cpp and inimate using haskell code
// Min Zhang
// Feb 11, 2016

#include <iostream>

using namespace std;

int main(){
  int age;
  int numberOfPeople = 0;
  int sumOfAge = 0;

  cout << "Input age, type -1 if you want to end the program." << endl;
  cin >> age;

  while (age != -1){
    sumOfAge += age;
    numberOfPeople ++;
    cout << "type another age, -1 to end: " << endl;
    cin >> age; 
  }

  cout << "total: " << numberOfPeople << " and average age is: " << sumOfAge/numberOfPeople << endl;
  
  return 0;
}
