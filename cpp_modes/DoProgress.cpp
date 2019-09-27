#include <iostream>
#include "wavefrf.h"
using namespace std;
void DoProgress(int step, int total, int& posplot){
	int posVec[12] = {1, 10 , 20, 30, 40, 50, 60, 70, 80, 90, 98, 100};
	int percent = (step * 100) / total;

	if (percent >= posVec[posplot] && posplot < 12) {
		cout << "We are about " << posVec[posplot] << "% done!" << endl;
		posplot++;
	}
}

