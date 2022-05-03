#pragma once
#include <iostream>
#include <vector>
#include <bitset> 

using namespace std;

class Subterranean {
public:
	/* Parameters in thetaOffset are the Theta shifting offsets. */
	vector<unsigned int> thetaOffset;

	/* Parameter piOffset is the Pi offset. */
	unsigned int piOffset = 12;

	/* Parameter piOffsetInverse is the inverse of Pi offset. */
	unsigned int piOffsetInverse = 150;

	/* Parameter Size represents the length of the state. */
	static const unsigned int Size = 257;

	/* Parameter thetaOffsetInverse represents the offsets of the inverse of Theta. */
	vector<unsigned int> thetaOffsetInverse{ 2, 4, 5, 8, 9, 12, 14, 16, 19, 20, 21, 22, 25, 26, 28, 29, 31, 37, 39, 42, 44,
 45, 49, 52, 53, 54, 58, 59, 60, 61, 62, 63, 64, 65, 71, 72, 73, 76, 77, 78, 79, 80, 82, 83, 86,
 90, 94, 95, 98, 99, 100, 102, 104, 105, 106, 108, 109, 111, 112, 119, 120, 124, 125, 127, 128,
 129, 130, 134, 136, 137, 138, 139, 141, 143, 145, 147, 148, 149, 150, 151, 152, 154, 158, 160,
 162, 163, 165, 166, 167, 172, 173, 174, 175, 177, 178, 179, 181, 184, 185, 187, 190, 193, 201,
 206, 209, 211, 216, 217, 219, 221, 222, 225, 226, 229, 231, 233, 236, 237, 238, 239, 242, 243, 245, 246, 248, 254, 256 };

	/** Parameter BeforeTheta shows the value of the state before Theta (after Chi). */
	bitset<Size> BeforeTheta;

	/** Vector AfterTheta shows the value of the state after Theta. */
	bitset<Size> AfterTheta;

	/** Parameter BeforeChi shows the value of the state before Chi (after Pi). */
	bitset<Size> BeforeChi;

public:
	/* Method to get all required parameters of Subterranean. */
	void GetParameters();

	/* Method to compute the inverse of the pi offset. */
	int InversePiParam();

	/* Method to compute the inverse of the Theta offsets. */
	void InverseThetaParam();

	/* Method to compute the value of the state after Theta given a state before Theta. */
	void Theta(bitset<Size> &arr, int aP);

	/* Method to compute the value of the state after Pi given a state before Pi. */
	void Pi(bitset<Size> &arr, int aP);

	/* Method to compute the value of the state before Theta given a state after Theta. */
	void ThetaInverse(bitset<Size> &input, int aP);

	/* Method to compute the value of the state before Pi given a state after Pi. */
	void PiInverse(bitset<Size> & input, int aP);
};