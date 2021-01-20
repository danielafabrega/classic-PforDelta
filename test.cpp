
//============================================================================
// Name        : TestPforDelta.cpp
// Author      : Hector Ferrada
// Copyright   : Your copyright notice
// Description : This use the PforDelta encoding method to compress an input array A[1..n]
//				 If you want to execute an exhaustive test, set PforDelta::TEST = true
//============================================================================

#include <string>
#include <iostream>
#include <sys/resource.h>
#include <sys/time.h>
#include <random>
#include "pfordelta.h"

using namespace std;

uint REPET = 1000000;

#define P 90 			// P% of probability that any number A[i] be in the range [L..R]
#define L 40
#define R 50
#define M 500		// Maximum value for any number of A

int main(int argc, char** argv)
{
	
	ulong i, n, prob;
	n = 40;
	ulong *A = new ulong[n];

	for (i=0; i<n ; i++){
		prob = rand()%100;
		if (prob < P)
			A[i] = L + rand()%(R-L);
		else
			A[i] = rand()%M;
	}
	A[2]=40;
	A[3]=64;
	A[4]=129;

	if (true){
		for (i=0; i<n; i++)
			cout << A[i] << " ";
		cout << endl;
	}
	

	cout << "Compressing array A[1.." << n << "] with PforDelta encoding... " << endl;
	cout<<"hola"<<endl;
	PforDelta *PfD = new PforDelta(A, n, 64);

	char directory[] = "./";

	//PfD->saveStructure(directory);
	PfD->testPforDelta(A);

	//PfD->~PforDelta();

	//PfD = new PforDelta(directory);

	cout << "################## " << endl;
	return 0;
}
