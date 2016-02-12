#include <iostream>
#include <cmath> // need power function
using namespace std;

// main1 print 1 to 9
int main1(){
  
  // iniitialization; loop end; increment
  for (int x = 1; x < 10; x++){
    cout << x << endl;
  }
  return 0;

}

// interest rate calculation
int main2(){

  float amount;
  float principle = 10000;
  float rate = 0.05;   

  for (int day = 1; day <= 30; day++){
    amount = principle * pow(1+rate, day);
    cout << day << " --- " << amount << endl;
  }
  return 0;
}

// do while loop; the purpose of having do-while loop
// the first number guarantee to run.
int main3(){
  int x;
  
  cin >> x;
  do{
  cout << x << endl;
  x++;
  }while(x<50);
  return 0;
}

// switch statement; replace if else statement
// BTW, does R have a switch?
int main(){
  int age;

  cin >> age;
  
  switch(age){
    case 16: 
      cout << "print 16" << endl;
      break; // if no break, will run default
    case 18:
      cout << "go buy lottery tickets" << endl;
      break;
    case 21:
      cout << "buy beer" << endl;
      break;
    default: 
      cout << "this is a default" << endl;
      // end of switch without break
  }
  return 0;
}
