/*classic PForDelta program, 
based on the algorithm of Dr. Hector Ferrada 
created on June 13, 2017
and PforDelta from chapter 9: Integer Encoding*/

#include <string>
#include <iostream>
#include <sys/resource.h>
#include <sys/time.h>
#include <random>
#include "pfordelta.h"

using namespace std;

/*PforDelta::PforDelta(char *pathFile){
	loadStructure(pathFile);
}
*/
PforDelta::PforDelta(ulong *A, ulong len, uint bitPerCell) {
	ulong i, j, k, l, m, max, num, auxMax;
	ulong bestBits, cont, cont2, totBits, diffLim, diffBit, diffEx;
    ulong limit, acum, blimit;
	uint lgNum, lgX, lgEMax;
	ulong *C, *Min, *Max;
	bool isEncoding = false;
	origBitA = bitPerCell;

	n = nP = len;
	nEMin = nEx = 0;
    auxMin, acum = 0;
	Ex = P = nullptr;
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
	/*
    cout<<"IMPRIMIR C: "<<endl;
	for (i=0; i<=maxBitsA; i++){
		cout<<"i: "<<i<<endl;
		cout<<C[i]<<endl;
		
	}*/
    //find 90% 
    limit = n*0.9;
    l=0;
    while(l<maxBitsA && acum<limit){
        //cout<<"C[l]"<<C[l]<<endl;
        acum +=C[l];
        l++;
    }
	//cout<<"l: "<<l<<endl;
	blimit = l-1;

	for (j=l; C[j]==0; j++);
	//cout<<"j: "<<j<<endl;
	if (C[j]){
		auxMin = Min[j];
	}
	//cout<<"auxMin: "<<auxMin<<endl;


   // cout<<"acum: "<<acum<<endl;
    
	//cout<<"basmin: "<<basMin<<endl;
	diffEx = limP - auxMin;
	diffLim = Max[blimit] - basMin;
	diffBit = 1 + (uint)(log(diffLim)/log(2));
	diffExBit = 1 + (uint)(log(diffEx)/log(2));
	cout<<diffExBit<<"diff"<<endl;
	//cout<<"diffBit"<<diffBit<<endl;
	b = diffBit;
	nP = acum;
	nEx = n-acum; //cantidad de números en las excepciones

	//ESTRUCTURA de P

	ulong bytesA = n*maxBitsA/8;
	BS = bit_vector(n, 0);
	
	k = nP*b / (8*sizeof(ulong));
	//cout<<"k: "<<k<<endl;
	if ((nP*b) % (8*sizeof(ulong))){
		k++;
	}
	//cout<<"k 2: "<<k<<endl;
	P = new ulong[k];
	sizePFD = k*sizeof(ulong);

	//Estructura de las excepciones
	if (nEx){
		BEx = bit_vector(n, 0);
		for (int m=0; m<n; m++){
		cout<<"Bex: "<<BEx[m]<<endl;
	}
		//cout<<"size; "<<BEx.size()<<endl;
		//cout<<"ENTRA"<<nEx<<endl;
		k = nEx*diffExBit / (8*sizeof(ulong));
		if ((nEx*diffExBit) % (8*sizeof(ulong))){
			k++;
		}
		cout<<"KKKKKK: "<<k<<endl;
		Ex = new ulong[k];
		sizePFD += k*sizeof(ulong);
		//cout<<"size: "<<sizePFD<<endl;
	}

	//ver si es exepcion o no excepcion
	ulong c = 0;
	ulong cEX = 0;
	for (i=j=k=0; i<n; i++, j+=bitPerCell){
		num = getNum64(A, j, bitPerCell);
		if (num <auxMin){
			BS[i] = 1;
			setNum64(P, k, b, num-basP);
			k += b;
		}
		else{
			//BEx[cEX] = 1;
			setNum64(Ex, c, diffExBit, num-auxMin);
			c += diffExBit;
			cout<<"c!!!!!!!!!!!: "<<c<<endl;
			cout<<"auxMIn!!!!!!!!!!!: "<<auxMin<<endl;
		}
		//cEX++;
		

	}
	for (int m=0; m<2;m++){
		cout<<"Ex :"<<Ex[m]<<endl;
	}
	for (int m=0; m<n; m++){
		cout<<"Bex: "<<BEx[m]<<endl;

	}
	BS_rrr = rrr_vector<127>(BS);
	sizePFD += size_in_bytes(BS_rrr);
	BS_ra = rrr_vector<127>::rank_1_type(&BS_rrr);

	if (nEx){
		BEx_ra = rank_support_v<>(&BEx);
		sizePFD += size_in_bytes(BEx);
	}

	decltype(BS) empty;
	BS.swap(empty);
    
}

void PforDelta::testPforDelta(ulong *A){
	ulong i, j, k, num, pos, *X;

	cout << "Testing Extract A[i] ..." << endl;
	for (i=0; i<n ; i++){
		num = extract(i);
		k = getNum64(A, i*origBitA, origBitA);
		cout<<"num: "<<num<<endl;
		if (num != k){ //compara si efectivamente los numeros que se extrayeron corresponden a los reales
			cout << "ERR. x = " << num << " != A[" << i << "] = " << k << endl;
			exit(1);
		}
		else{
			cout<<"OKKKKK"<<endl;
		}
	}/*
	uint bitX, inc = 10;
	cout << "Testing Extract A[i..i+" << inc << "] ..." << endl;
	if (n > 10*inc){
		cout<<"NO DEBERIA ENTRAR"<<endl;
		for (i=0; i<n-inc ; i+=inc/2){
			extract(i, i+inc, &X, &bitX);
			for (k=i, pos=0; k<i+inc; k++, pos+=bitX){
				num = getNum64(X, pos, bitX);
				cout<<"num: "<<num<<endl;
				j = getNum64(A, k*origBitA, origBitA);
				if (num != j){
					cout << "ERR. x = " << num << " != A[" << k << "] = " << j << endl;
					exit(1);
				}
			}
		}
	}else{
		cout<<"AQUI SIII"<<endl;
		extract(0, n-1, &X, &bitX);
		for (k=0, pos=0; k<n; k++, pos+=bitX){
			num = getNum64(X, pos, bitX);
			cout<<"num: "<<num<<endl;
			j = getNum64(A, k*origBitA, origBitA);
			if (num != j){
				cout << "ERR. x = " << num << " != A[" << k << "] = " << j << endl;
				exit(1);
			}
		}
	}*/
	cout << "Test OK !!" << endl;
}
ulong PforDelta::getNum64(ulong *A, ulong ini, uint len){
	ulong i=ini>>BW64, j=ini-(i<<BW64);
	ulong result = (A[i] << j) >> (W64-len);
	if (j+len > W64)
		result = result | (A[i+1] >> (WW64-j-len));

	return result;
}

// set the number x as a bitstring sequence in *A. In the range of bits [ini, .. ini+len-1] of *A. Here x has len bits
void PforDelta::setNum64(ulong *A, ulong ini, uint len, ulong x) {
	ulong i=ini>>BW64, j=ini-(i<<BW64);
	//cout<<"i: "<<i<<" j: "<<j<<endl;
	//cout<<"AANTES A[ "<<i<<"]: "<<A[i]<<endl;
	//ulong myMask;
	if ((j+len)>W64){
		//cout<<"J: "<<j<<endl;
		//cout<<"ENTRA W64"<<endl;
		//cout<<"myMask: "<<myMask<<endl;
		ulong myMask = ~(~0ul >> j); //~0ul = 2 elevado a 64, deja solo los bits superiores 
		//cout<<"myMask: "<<myMask<<endl;
		
		//cout<<"AANTES A[i]: "<<A[i]<<endl;
		A[i] = (A[i] & myMask) | (x >> (j+len-W64));
		//cout<<"A[i]: "<<A[i]<<endl;
		myMask = ~0ul >> (j+len-W64);
		A[i+1] = (A[i+1] & myMask) | (x << (WW64-j-len));
	}else{
		//cout<<"ENTRA W64 ELSE"<<endl;
		//cout<<"J else: "<<j<<endl;
		ulong myMask = (~0ul >> j) ^ (~0ul << (W64-j-len)); // XOR: 1^1=0^0=0; 0^1=1^0=1
		//cout<<"myMask else: "<<myMask<<endl;
		A[i] = (A[i] & myMask) | (x << (W64-j-len));
	}
}

ulong PforDelta::extract(ulong i){
	ulong r = BS_ra.rank(i+1);
	cout<<"r: "<<r<<endl;

	if (BS_rrr[i])
		return basP+getNum64(P, (r-1)*b, b);
	else{
		//cout<<"ERROR"<<endl;
		
		ulong ze = i+1-r;
		//cout<<"ERROR 222"<<endl;
		return auxMin+getNum64(Ex, (ze-1)*diffExBit, diffExBit);
		//cout<<"GET NUMBER!!!!!!!!!!!!!!!!!"<<getNum64(Ex, 0, 8)<<endl;
		//return auxMin+getNum64(Ex, 0, 8);
		
	return 0;
	}
}

/*void PforDelta::extract(ulong i, ulong j, ulong **X, uint *nBit){
	ulong x, k, one, oneEx;

	k = maxBitsA*(j-i+1) / (8*sizeof(ulong));
	if ((maxBitsA*(j-i+1)) % (8*sizeof(ulong)))
		k++;
	//cout<<"KKK: "<<k<<endl;
	ulong *AUX = new ulong[k];

	one = k = 0;
	if (i){
		one = BS_ra.rank(i);	// ones previous to i
		k = i-one;				// number of zeros previous to i
	}

	if (nEx){
		if (k)
			oneEx = BEx_ra.rank(k);
		else
			oneEx = 0;

		one *= b;
		for (ulong pos=0; i<=j; i++, pos+=maxBitsA){
			if (BS_rrr[i]){
				//cout<<"ENTRA BS_RRR"<<endl;
				x = basP+getNum64(P, one, b);
				one+=b;
			}else{
				if (BEx[k]){
					x = auxMin+getNum64(Ex, (k-oneEx)*diffExBit, diffExBit);
				}
				k++;
			}
			//cout<<"XXXX = "<<x<<endl;
			setNum64(AUX, pos, maxBitsA, x);
		}

	*nBit = maxBitsA;
	*X = AUX;
	/*for (int i=0; i<=k; i++){
		cout<<AUX[i]<<"AUXA"<<endl;
	}
	}
}*/

/*
void PforDelta::saveStructure(char *pathName){
	char *fileName = new char[400];

	strcpy(fileName, "");
	strcpy(fileName, pathName);
	strcat(fileName, "PforD_ST.pfd");
	ofstream os (fileName, ios::binary);
	cout << "Save PforDelta data structures in " << fileName << endl;

	os.write((const char*)&n, sizeof(ulong));
	os.write((const char*)&nP, sizeof(ulong));
	os.write((const char*)&nEx, sizeof(ulong));

	os.write((const char*)&basMin, sizeof(ulong));
	os.write((const char*)&basP, sizeof(ulong));
	os.write((const char*)&limP, sizeof(ulong));

	os.write((const char*)&origBitA, sizeof(uint));
	os.write((const char*)&minBitsA, sizeof(uint));
	os.write((const char*)&maxBitsA, sizeof(uint));
	os.write((const char*)&b, sizeof(uint));
	os.write((const char*)&diffExBit, sizeof(uint));


	ulong lenArray = nP*b / (8*sizeof(ulong));
	if ((nP*b) % (8*sizeof(ulong)))
		lenArray++;
	os.write((const char*)P, lenArray*sizeof(ulong));
	cout << " .- P[] " << lenArray*sizeof(ulong) << " Bytes" << endl;


	if (nEx){
		lenArray = nEx*diffExBit / (8*sizeof(ulong));
		if ((nEx*diffExBit) % (8*sizeof(ulong)))
			lenArray++;
		os.write((const char*)Ex, lenArray*sizeof(ulong));
		cout << " .- ExMax[] " << lenArray*sizeof(ulong) << " Bytes" << endl;
	}

	os.close();

	strcpy(fileName, "");
	strcpy(fileName, pathName);
	strcat(fileName, "PforD_BS_rrr.pfd");
	store_to_file(BS_rrr, fileName);
	cout << " .- BS_rrr " << size_in_bytes(BS_rrr) << " Bytes" << endl;

	strcpy(fileName, "");
	strcpy(fileName, pathName);
	strcat(fileName, "PforD_BS_ra.pfd");
	store_to_file(BS_ra, fileName);

	if(nEx){
		strcpy(fileName, "");
		strcpy(fileName, pathName);
		strcat(fileName, "PforD_BEx.pfd");
		store_to_file(BEx, fileName);

		strcpy(fileName, "");
		strcpy(fileName, pathName);
		strcat(fileName, "PforD_BEx_ra.pfd");
		store_to_file(BEx_ra, fileName);
	}

	cout << "   PforDelta Structures saved !" << endl;
	cout << "______________________________________________________________" << endl;
}

void PforDelta::loadStructure(char *pathName){
	char *fileName = new char[400];
	ulong lenArray;

	strcpy(fileName, "");
	strcpy(fileName, pathName);
	strcat(fileName, "PforD_ST.pfd");
	ifstream is(fileName, ios::binary);
	cout << " Load data structure from " << fileName << endl;

	is.read((char*)&n, sizeof(ulong));
	is.read((char*)&nP, sizeof(ulong));
	is.read((char*)&nEx, sizeof(ulong));

	is.read((char*)&basMin, sizeof(ulong));
	is.read((char*)&basP, sizeof(ulong));
	is.read((char*)&limP, sizeof(ulong));

	is.read((char*)&origBitA, sizeof(uint));
	is.read((char*)&minBitsA, sizeof(uint));
	is.read((char*)&maxBitsA, sizeof(uint));
	is.read((char*)&b, sizeof(uint));
	is.read((char*)&diffExBit, sizeof(uint));


	lenArray = nP*b / (8*sizeof(ulong));
	if ((nP*b) % (8*sizeof(ulong)))
		lenArray++;
	P = new ulong[lenArray];
	is.read((char*)P, lenArray*sizeof(ulong));
	cout << " ** size of P[] " << lenArray*sizeof(ulong) << " Bytes" << endl;
	sizePFD = lenArray*sizeof(ulong);


	if (nEx){
		lenArray = nEx*diffExBit/ (8*sizeof(ulong));
		if ((nEx*diffExBit) % (8*sizeof(ulong)))
			lenArray++;
		Ex = new ulong[lenArray];
		is.read((char*)Ex, lenArray*sizeof(ulong));
		cout << " ** size of ExMax[] " << lenArray*sizeof(ulong) << " Bytes" << endl;
		sizePFD += lenArray*sizeof(ulong);
	}
	is.close();

	strcpy(fileName, "");
	strcpy(fileName, pathName);
	strcat(fileName, "PforD_BS_rrr.pfd");
	load_from_file(BS_rrr, fileName);
	cout << " ** size of BS_rrr " << size_in_bytes(BS_rrr) << " Bytes" << endl;
	sizePFD += size_in_bytes(BS_rrr);

	strcpy(fileName, "");
	strcpy(fileName, pathName);
	strcat(fileName, "PforD_BS_ra.pfd");
	load_from_file(BS_ra, fileName);
	util::init_support(BS_ra, &BS_rrr);

	if(nEx){
		strcpy(fileName, "");
		strcpy(fileName, pathName);
		strcat(fileName, "PforD_BEx.pfd");
		load_from_file(BEx, fileName);
		sizePFD += size_in_bytes(BEx);

		strcpy(fileName, "");
		strcpy(fileName, pathName);
		strcat(fileName, "PforD_BEx_ra.pfd");
		load_from_file(BEx_ra, fileName);
		sizePFD += size_in_bytes(BEx_ra);
		util::init_support(BEx_ra, &BEx);
	}

	cout << "   PforDelta Structures loaded. sizePFD = " << sizePFD << " Bytes." << endl;
	cout << "______________________________________________________________" << endl;

}

// ====================================================================================================
PforDelta::~PforDelta() {
	delete [] P;


	if (nEx){
		delete [] Ex;
		cout<<"Erro1"<<endl;
	}

	decltype(BS_rrr) empty;
	BS_rrr.swap(empty);
	cout<<"Erro2"<<endl;

	if (nEx){
		cout<<"Erro3"<<endl;

		decltype(BEx) emptyBEx;
		BEx.swap(emptyBEx);
		cout<<"Erro4"<<endl;

		decltype(BEx_ra) emptyBEx_ra;
		BEx_ra.swap(emptyBEx_ra);
		cout<<"Erro5"<<endl;

	}
}
*/