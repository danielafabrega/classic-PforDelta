/*
 * PforDelta.h
 *
 *  Created on: 13 Jan 2017
 *      Author: Hector Ferrada
 */

#ifndef PFORDELTA_H_
#define PFORDELTA_H_

//#include <sdsl/bit_vectors.hpp>

//using namespace sdsl;

using namespace std;

typedef unsigned int uint;
typedef unsigned long ulong;


const uint W64 = 64;
const uint BW64 = 6;	// pow of two for W64
const uint WW64 = 128;

class PforDelta {
public:
	static bool TRACE;	// true: print all details for console
	static bool TEST;	// true: print all details for console

	//bit_vector BS;							// BitVector of length n to mark numbers in the PforDelta interval
	//rrr_vector<127> BS_rrr;
	//rrr_vector<127>::rank_1_type BS_ra;

	//bit_vector BEx;							// BitVector of length nP to mark exceptions
	//rank_support_v<> BEx_ra;

	ulong basMin;
	ulong *Ex;
	ulong basP;
	ulong *P;
	ulong limP;
	ulong *ExMax;
	ulong n;		// = nS + nM + nL
	ulong nEMin;	// # 1's in BEx
	ulong nP;		// # 1's in BS
	ulong nEx;	// # 0's in BEx

	uint origBitA; 	// bits for each value in the original A[] array received in the parameter 'bitPerCell'
	
	uint minBitsA;	// bits for the minimum value in the original A[] array
	uint maxBitsA;	// bits for the maximum value in the original A[] array
	uint bEMin;		// bits for each small number in ExL[]
	uint b;			// bits for each interval item in P[]
	uint bEMax;		// bits for each large number in ExR[]

	ulong sizePFD;
	bool isBEx;			// it is true if we create the bit_vector BEx and its rank support

	PforDelta(ulong *A, ulong n, uint bitPerCell);
	//PforDelta(char *pathFile);
	//void testPforDelta(ulong *A);

	// retrieve A[i]
	//ulong extract(ulong i);
	// retrieve A[i..j]
	//void extract(ulong i, ulong j, ulong **X, uint *nBit);

	// set the number x as a bitstring sequence in *A. In the range of bits [ini, .. ini+len-1] of *A. Here x has len bits
	//void setNum64(ulong *A, ulong ini, uint len, ulong x);

	// return (in a unsigned long integer) the number in A from bits of position 'ini' to 'ini+len-1'
	ulong getNum64(ulong *A, ulong ini, uint len);

	// save load structure to/from pathFile
	//void saveStructure(char *pathFile);
	//void loadStructure(char *pathFile);

	//virtual ~PforDelta();
};


#endif /* PFORDELTA_H_ */
