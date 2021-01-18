/*classic PForDelta program, 
based on the algorithm of Dr. Hector Ferrada 
created on June 13, 2017*/

#include <string>
#include <iostream>
#include <sys/resource.h>
#include <sys/time.h>
#include <random>
#include "pfordelta.h"

using namespace std;


PforDelta::PforDelta(ulong *A, ulong len, uint bitPerCell) {
	ulong i, j, k, l, m, max, num, auxMin, auxMax;
	ulong bestBits, cont, cont2, totBits, diffLim, diffBit;
    ulong limit, acum, blimit;
	uint lgNum, lgX, lgEMax;
	ulong *C, *Min, *Max;
	bool isEncoding = false;
	isBEx = false;
	origBitA = bitPerCell;

	n = nP = len;
	nEMin = nEx = 0;
    acum = 0;
	ExMin = ExMax = P = nullptr;
	basMin = max = getNum64(A, 0, origBitA);
	for (i=1, j=origBitA; i<n; i++, j+=origBitA){
		num = getNum64(A, j, origBitA);
		
		if (num < basMin)
			basMin = num;
		if (num > max)
			max = num;
	}
	minBitsA = 1 + (uint)(log(basMin)/log(2));
	maxBitsA = 1 + (uint)(log(max)/log(2));
	basP = basMin; 
	limP = max; 
	b = maxBitsA;

	C = new ulong[maxBitsA+1];
	Min = new ulong[maxBitsA+1];
	Max = new ulong[maxBitsA+1];
	for (i=0; i<=maxBitsA; i++){
		C[i] = Max[i] = 0;
		Min[i] = max;
	}

	for (i=j=0; i<n; i++, j+=origBitA){
		num = getNum64(A, j, origBitA);
		//cout<<"num: "<<num<<endl;
		lgNum = 1 + (uint)(log(num)/log(2));
		(C[lgNum])++;

		if (num < Min[lgNum])
			 Min[lgNum] = num;
		if (num > Max[lgNum])
			 Max[lgNum] = num;
	}
    cout<<"IMPRIMIR C: "<<endl;
	for (i=0; i<=maxBitsA; i++){
		cout<<"i: "<<i<<endl;
		cout<<C[i]<<endl;
		
	}
    //find 90% 
    limit = n*0.9;
    l=0;
    while(l<maxBitsA && acum<limit){
        //cout<<"C[l]"<<C[l]<<endl;
        acum +=C[l];
        l++;
    }
	cout<<"l: "<<l<<endl;
	//preguntar si el siguiente 
	
	// avanzar hasta que sea dif de 0

    cout<<"acum: "<<acum<<endl;
    blimit = l-1;
	cout<<"basmin: "<<basMin<<endl;
	diffLim = Max[blimit] - basMin;
	diffBit = 1 + (uint)(log(diffLim)/log(2));
	cout<<"diffBit"<<diffBit<<endl;
	b = diffBit;
	nP = acum;
	nEx = n-acum; //cantidad de números en las excepciones

	//ESTRUCTURA

	ulong bytesA = n*maxBitsA/8;
	k = nP*b / (8*sizeof(ulong));
	cout<<"k: "<<k<<endl;
	if ((nP*b) % (8*sizeof(ulong))){
		k++;
	}
	cout<<"k 2: "<<k<<endl;
	P = new ulong[k];
	sizePFD = k*sizeof(ulong);

	/*if (nEx){
		cout<<"ENTRA A nEMAX: "<<endl;
		k = nEx*bEMax / (8*sizeof(ulong)); 
		if ((nEx*bEMax) % (8*sizeof(ulong)))
			k++;
		ExMax = new ulong[k];
		sizePFD += k*sizeof(ulong);
	}*/
    
}
ulong PforDelta::getNum64(ulong *A, ulong ini, uint len){
	ulong i=ini>>BW64, j=ini-(i<<BW64);
	ulong result = (A[i] << j) >> (W64-len);
	if (j+len > W64)
		result = result | (A[i+1] >> (WW64-j-len));

	return result;
}

