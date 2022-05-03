#include "Subterranean.h"
#include <iostream>
#include <stdio.h>
#include <vector>
#include <algorithm>


using namespace std;

void Subterranean::GetParameters() {
	/* The offset values of the Theta from the smallest one to the largest one. */ 
	thetaOffset.push_back(0);
	thetaOffset.push_back(3);
	thetaOffset.push_back(8);
	/* The offset value of the Pi. */
	piOffset = 12;
}

int Subterranean::InversePiParam() {
	for (piOffsetInverse = 2; piOffsetInverse < Size; piOffsetInverse++) {
		/* If the offset value of Pi function is equal to its inverse. */
		if ((piOffset * piOffsetInverse) % Size == 0) {
			cout << "The pi parameter is not well chosen!\n"
				<< "GCD(" << piOffset << "," << Size << ") = 0 is not acceptable.\n";
		}
		/* Returns the inverse of Pi offset. */
		else if ((piOffset * piOffsetInverse) % Size == 1) {
			return piOffsetInverse;
		}
	}
	return 0;
}

void Subterranean::InverseThetaParam() {
	//int inv[Size + 1] = { 0 };
	//inv[Size] = 1;
	//for (unsigned int i = Size; i >= 0; i--) {
	//	if (inv[i] == 1) {
	//		thetaOffsetInverse.push_back(i - thetaOffset.back());//245
	//		for (unsigned int j = 0; j < thetaOffset.size(); j++) {
	//			inv[thetaOffsetInverse.back() + thetaOffset[j]] ^= 1;//253,251
	//		}
	//	}
	//}
	//sort(thetaOffsetInverse.begin(), thetaOffsetInverse.end());
	;
}

void Subterranean::Theta(bitset<Size> &arr, int aP) {
	for (unsigned int i = 0; i < thetaOffset.size(); i++) {
		arr.flip((aP - thetaOffset[i] + Size) % Size);
	}
}

void Subterranean::Pi(bitset<Size>& arr, int aP) {
	for (unsigned int i = 0; i < thetaOffset.size(); i++) {
		arr.flip((((aP - thetaOffset[i] + Size) % Size) * piOffsetInverse) % Size);
	}
}

void Subterranean::PiInverse(bitset<Size>& input, int aP) {
	input.flip((aP * piOffset) % Size);
}

void Subterranean::ThetaInverse(bitset<Size>& input, int aP) {
	for (unsigned int j = 0; j < thetaOffsetInverse.size(); j++) {
		input.flip(((aP * piOffset) + thetaOffsetInverse[j] + Size) % Size);
	}
}


