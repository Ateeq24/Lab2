#include"Mat.h"
#include <ctime>
#include <iostream>

using namespace std;

int main()
{
	Test test1;
	if (test1.initialize_test())
		cout << "Successfully initialized Matrix" << endl;
	else
		cout << "Matrix initialization Failed!" << endl;

	time_t start = time(NULL);
	start *= 1000;

	if (test1.Reg_Mul_Test())
		cout << "Regular Multiplication Test Passed!" << endl;
	else
		cout << "Regular Multiplication Test Failed!" << endl;

	time_t final = time(NULL);
	final *= 1000;

	cout << "Total Time: " << final - start << endl;


	start = time(NULL);
	if (test1.Strassen_Mul_Test())
		cout << "Multiplication through Stressen Algorithm Test Passed!" << endl;
	else
		cout << "Multiplication through Stressen Algorithm Test Failed!" << endl;

	final = time(NULL);

	cout << "Total Time: " << final - start << endl;
}

