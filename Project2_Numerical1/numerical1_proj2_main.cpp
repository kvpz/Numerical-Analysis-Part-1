#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
using namespace std;

void nfactorial(int n, double& nfact) 
{
	int i;
	nfact = 1.0;

	for (i = 1; i <= n; i++) 
	{
		nfact = nfact*i;
	}
}
/*
double power(double x, int N) 
{
	double pow = 1.0;
	for (int i = 1; i <= N; i++)
		pow *= x;
	cout << "current power value: " << pow << endl;
	return pow;
}
*/

void singlePrecision()
{
	float x = 1.0f;
	float h = 0.1f;//reset value of h
	float htemp = 0.1f;
	float M = 0.0f; //reset value of M

	cout << setprecision(20);
	cout << "Direct implementation:" << endl;
	//for exact values of M(h), M(0) formula
	for (int n = 1; n <= 20; n++)
	{
		M = (sqrt(x + h) - sqrt(x)) / h; //calculates M(h) directly
		cout << "M(" << h << ") = " << M << endl;
		h *= htemp;
	}//for loop



}

void doublePrecision()
{
	double x = 1.0;
	double h = .1;//reset value of h
	double htemp = 0.1;
	double M = 0.0;

	cout << setprecision(20);
	cout << "Direct implementation:" << endl;
	//for exact values of M(h), M(0) formula
	for (int n = 1; n <= 20; n++)
	{
		M = (sqrt(x + h) - sqrt(x)) / h; //calculates M(h) directly
		cout << "M(" << h << ") = " << M << endl;
		h *= htemp;
	}//for loop

}



//output to file
void outputLogSp(ostream& fileSp)
{ //assume x = 1
	float x = 1.0f;
	float h = 0.1f;
	float htemp = 0.1f;
	float Mzero = 0.5f;
	float Mh = 0.0f;

	for (int i = 1; i < 21; i++)
	{
		Mh = 1/(2*sqrt(x+h));
		cout << "MH: " << Mh << endl;
		fileSp << log(h) << ' ' << log(fabs(Mh - Mzero)) << endl;
		h *= htemp; //h*=0.1 is not allowed
	}


}

void outputLogDp(ostream& fileDp)
{
	double x = 1.0;
	double h = 0.1;
	double htemp = 0.1;
	double Mzero = 0.5;
	double Mh = 0.0;
	cout << setprecision(15);
	for (int i = 1; i < 21; i++)
	{
		Mh = 1/(2*sqrt(x+h));
		cout << "MH: " << Mh << endl;
		fileDp << log(h) << ' ' << log(fabs(Mh - Mzero)) << endl;
		h *= htemp; //h*=0.1 is not allowed
	}

}

void menu()
{
	cout << "To output M(h) for each h^-x for x=1,...,20 (as h approaches 0)" << endl;
	cout << "Type in number to select an option:" << endl;
	cout << "1) singlePrecision" << endl;
	cout << "2) doublePrecision" << endl;
	cout << "3) output to file (single precision)" << endl;
	cout << "4) output to file (double precision)" << endl;
	cout << "0) End program" << endl;
}

void selection(ofstream& fileSp,ofstream& fileDp, bool& end)
{
	int	c = 0;
	cin >> c;
	switch (c)
	{
	case 0:
	{ //"non-labels" must be in brackets
		end = true;
		break;
	}
	case 1:
		singlePrecision();
		break;
	case 2:
		doublePrecision(); 
		break;
	case 3:
		outputLogSp(fileSp);
		break;
	case 4:
		outputLogDp(fileDp);
		break;
	default:
		menu();
	}
}

int main()
{
	
	ofstream fileSp;
	ofstream fileDp;
	fileSp.open("spOutput_Num1_proj2.txt");
	fileDp.open("dpOutput_Num1_Proj2.txt");
	bool end = false; //if true, program terminates

	float a = 0.0f;
	float b = 0.5f;
	cout << a - b << endl;

	do{
		menu();
		selection(fileSp,fileDp, end);

	} while (end == false);

	fileSp.close();
	fileDp.close();
	return 0;
}

/*
NOTES
floating point values stored in memory are much more precise because there are
many more floating point registers. This is especially when dealing with floating double precision.
http://stackoverflow.com/questions/7972634/single-precision-and-double-precision-casting-gives-different-answers

double has 15 to 16 decimal digits of precision, while float only has 7

h*=h IS NOT THE SAME AS h=h*.1 or h*=.1


This yields a less accurate result
temp = sqrt(x + h)/h;
temp = temp - sqrt(x)/h;

when you say 2nd order, must that also include the truncation error?

Is M(0), M(h) both only to be represented as taylor expansions? (is direct implementation necessary?)

c++ does not allow to send float as a template argument.

The limit as h->0 for (sqrt(x+h)-sqrt(x))/h = 1/(2*sqrt(x))
	x = 1, then the limit is 1/2=0.5.

	On a 32 bit machine, floats only have a precision of 6 to 7 digits and 
	values as high as 3.4028235E+38 and as low as -3.4028235E+38.
*/