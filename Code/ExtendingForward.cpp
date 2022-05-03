#pragma once
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdint.h>
#include <algorithm>
#include <vector> 
#include "Subterranean.h"
#include "TrailSeedIterator.h"

using namespace std;

void TrailSeedIterator::ForwardExtension() {
	setBasisOffset();

	if (ComputeLastRoundScore(NumberOfRounds) == true) {
		copyStateForward(NumberOfRounds);
		/* Generate affine space by adding basis vectors to the canonical offset. */
		/* 128 is the maximum number of bases or maximum weight of the last round. 
		* To be more precise, it should be equal to (BasisBeforeThetaF.size() / 3) instead of 128,
		* but, 128 is the maximum that guarantees we do not miss any output difference. */
		bitset<128> grayCounter = 0;
		for (int s = 0; s < (1 << ((BasisBeforeThetaF.size() / 3) - 1)); s++) {
			/*The way to generate the next value by gray code,
			* when we flip the bit at position "i + 1" if "i" is the LSB.*/
			for (unsigned int i = 0; i < grayCounter.size(); i++) {
				if (grayCounter[i] == 1) {
					grayCounter.flip(i + 1);
					updateOffset(i + 1);
					i = grayCounter.size() - 1;
				}
			}
			/*The way to generate the next value by gray code
			* when we flip the bit at position "0".*/
			grayCounter.flip(0);
			updateOffset(0);

		}
	}
}

void TrailSeedIterator::copyStateForward(int aNumberOfRounds) {
	if (aNumberOfRounds == 3) {
		if (ComputeMinRevWeight(TempBefTheta) + ComputeWeight(TempBefChi) + ComputeWeight(OffsetBeforeChiF) <= MaxCost4Round - Min_W2) {
			OffsetAfterThetaF.reset();
			for (int i = 0; i < Size; i++) {
				if (OffsetBeforeThetaF[i] == 1) {
					Theta(OffsetAfterThetaF, i);
				}
			}
			copyState(BeforeThetaR1, TempBefTheta);
			copyState(AfterThetaR1, TempAftTheta);
			copyState(BeforeChiR1, TempBefChi);
			copyState(BeforeThetaR2, OffsetBeforeThetaF);
			copyState(AfterThetaR2, OffsetAfterThetaF);
			copyState(BeforeChiR2, OffsetBeforeChiF);
			showCout();
		}
	}
	else {
		if (ComputeMinRevWeight(BeforeTheta) + ComputeWeight(BeforeChi) + ComputeWeight(TempBefChi) + ComputeWeight(OffsetBeforeChiF) <= MaxCost4Round) {
			copyState(BeforeThetaR3, OffsetBeforeThetaF);
			copyState(AfterThetaR3, OffsetAfterThetaF);
			copyState(BeforeChiR3, OffsetBeforeChiF);
			showCout4();

			copyState(BeforeThetaFinalR1, BeforeTheta);
			copyState(AfterThetaFinalR1, AfterTheta);
			copyState(BeforeChiFinalR1, BeforeChi);
			copyState(BeforeThetaFinalR2, TempBefTheta);
			copyState(AfterThetaFinalR2, TempAftTheta);
			copyState(BeforeChiFinalR2, TempBefChi);
			copyState(BeforeThetaFinalR3, OffsetBeforeThetaF);
			copyState(AfterThetaFinalR3, OffsetAfterThetaF);
			copyState(BeforeChiFinalR3, OffsetBeforeChiF);
		}
	}
}

void TrailSeedIterator::updateOffset(int aP) {
	/* There are three types of bases:
	* basis with one active bits and no overlap,
	* basis with two active bits and no overlap,
	* basis with two active bits and an overlap. */
	/* For all three types we have: */
	OffsetBeforeThetaF.flip(BasisBeforeThetaF[aP * 3 + 1]);
	for (unsigned int j = 0; j < thetaOffset.size(); j++) {
		OffsetBeforeChiF.flip(BasisBeforeChiF[aP * 2 * thetaOffset.size() + j]);
	}
	/* If it is not type1: */
	if (BasisBeforeThetaF[aP * 3] != 1) {
		OffsetBeforeThetaF.flip(BasisBeforeThetaF[aP * 3 + 2]);
		for (unsigned int k = thetaOffset.size(); k < 2 * thetaOffset.size(); k++) {
			OffsetBeforeChiF.flip(BasisBeforeChiF[aP * 2 * thetaOffset.size() + k]);
		}
	}
	copyStateForward(NumberOfRounds);
}

void TrailSeedIterator::setBasisOffset() {
	OffsetBeforeThetaF.reset();
	OffsetBeforeChiF.reset();
	//type, pos1, pos2
	BasisBeforeThetaF.clear();
	BasisAfterThetaF.clear(); 
	BasisBeforeChiF.clear();

	allOutputActiveBits.clear();

	for (int i = 0; i < Size; i++) {
		if (TempBefChi[i] == 1 && TempBefChi[(i + 1) % Size] == 0 && TempBefChi[(i + 2) % Size] == 0) {
			OffsetBeforeThetaF[i] = 1;
			Pi(OffsetBeforeChiF, i);
			allOutputActiveBits.push_back(i);
		}
		else if (TempBefChi[i] == 1 && TempBefChi[(i + 1) % Size] == 1 && TempBefChi[(i + 2) % Size] == 0 && TempBefChi[(i + 3) % Size] == 1) {
			/* We DO NOT add active bit to the allOutputActiveBits
			* instead, we will add them later in GenerateBases(). */
			OffsetBeforeThetaF[i] = 1;
			Pi(OffsetBeforeChiF, i);
		}

		if (TempBefChi[(i + 1) % Size] == 1 && TempBefChi[(i + 2) % Size] == 0 && TempBefChi[(i + 3) % Size] == 0) {
			GenerateBases(1, i);
		}
		else if (TempBefChi[(i + 1) % Size] == 1 && TempBefChi[(i + 2) % Size] == 1) {
			GenerateBases(1, i);
		}
		else if (TempBefChi[i] == 0 && TempBefChi[(i + 1) % Size] == 0 && TempBefChi[(i + 2) % Size] == 1) {
			GenerateBases(1, i);
		}
		else if (TempBefChi[(i + 1) % Size] == 1 && TempBefChi[(i + 2) % Size] == 0 && TempBefChi[(i + 3) % Size] == 1) {
			if (TempBefChi[i] == 0) {
				GenerateBases(2, i);
			}
			else {
				GenerateBases(3, i);
			}
		}
	}
}

void TrailSeedIterator::GenerateBases(int aType, int aP) {
	/* There are three types of bases:
	* type1: basis with one active bits. There is no overlap.
	* type2: basis with two active bits. There is no overlap.
	* type3: basis with two active bits. There is an overlap. */
	BasisBeforeThetaF.push_back(aType);//type
	BasisBeforeThetaF.push_back(aP);//pos1
	BasisBeforeThetaF.push_back((aP + 1) % Size);//pos2
	allOutputActiveBits.push_back(aP);
	if (aType != 1) {
		allOutputActiveBits.push_back((aP + 1) % Size);
	}
	for (int i = 0; i < 2; i++) {
		for (unsigned int m = 0; m < thetaOffset.size(); m++) {
			BasisAfterThetaF.push_back((aP + i - thetaOffset[m] + Size) % Size);
			BasisBeforeChiF.push_back(((aP + i - thetaOffset[m] + Size) * piOffsetInverse) % Size);
		}
	}
}

bool TrailSeedIterator::ComputeLastRoundScore(int aNumberOfRounds) {
	/* The score of the last round. */
	int LastRoundScore = 0;
	sort(allOutputActiveBits.begin(), allOutputActiveBits.end());
	//First and last active bits.
	if (allOutputActiveBits[0] - allOutputActiveBits.back() + Size > thetaOffset.back()) {
		//If the distance is greater than thetaOffset.back(), it adds at least 3 active bits.
		LastRoundScore += 3;
		if (allOutputActiveBits[0] - allOutputActiveBits.back() + Size > (2 * piOffset + thetaOffset.back())) {
			//If the distance is greater than (2 * piOffset + thetaOffset.back()), it adds at least 3 001-string.
			LastRoundScore += 3;
		}
	}

	//Other active bits.
	for (unsigned int i = 0; i < allOutputActiveBits.size() - 1; i++) {
		if (allOutputActiveBits[i + 1] - allOutputActiveBits[i] > thetaOffset.back()) {
			//If the distance is greater than thetaOffset.back(), it adds at least 3 active bits.
			LastRoundScore += 3;
			if (allOutputActiveBits[i + 1] - allOutputActiveBits[i] > (2 * piOffset + thetaOffset.back())) {
				//If the distance is greater than (2 * piOffset + thetaOffset.back()), it adds at least 3 001-string.
				LastRoundScore += 3;
			}
		}
	}

	if (aNumberOfRounds == 3) {
		if (ComputeMinRevWeight(TempBefTheta) + ComputeWeight(TempBefChi) + LastRoundScore + Min_W2 > MaxCost4Round) {
			return false;
		}
	}
	else {
		if (ComputeMinRevWeight(BeforeTheta) + ComputeWeight(BeforeChi) + ComputeWeight(TempBefChi) + LastRoundScore > MaxCost4Round) {
			return false;
		}
	}
	return true;
}