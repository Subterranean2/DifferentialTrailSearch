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

void TrailSeedIterator::BackwardExtension() {
	BeforeChiBack.reset();
	pos.clear();
	GenerateAllInputDiffInitialize();
	/* Generate the next input difference.*/
	while (BackwardExtensionNext() == true) {
		;
	}
}

void TrailSeedIterator::BackwardUpdateState() {
	if (ComputeWeight(BeforeChiBack) >= Min_W2 && ComputeWeight(BeforeChiBack) <= maxBudgetBackward && isCompatible(Size - 1, 0) == true) {
		BeforeThetaBack.reset();
		AfterThetaBack.reset();
		/* call pi and theta inverse to generate state BeforeThetaBack. */
		for (unsigned int i = 0; i < pos.size(); i++) {
			PiInverse(AfterThetaBack, pos[i]);
			ThetaInverse(BeforeThetaBack, pos[i]);
		}
		copyStateBackward(NumberOfRounds);
	}
}

void TrailSeedIterator::GenerateAllInputDiffInitialize() {
	BeforeThetaBack.reset();
	AfterThetaBack.reset();
	BeforeChiBack.reset();
	number = 1;
	pos.push_back(Size - 1);
	BeforeChiBack.flip(pos[number - 1]);
	BackwardUpdateState();
}

void TrailSeedIterator::copyStateBackward(int aNumberOfRounds) {
	if (aNumberOfRounds == 3) {
		if (ComputeMinRevWeight(BeforeThetaBack) + ComputeWeight(BeforeChiBack) + ComputeWeight(TempBefChi) <= MaxCost4Round - Min_W0) {
			copyState(BeforeThetaR1, BeforeThetaBack);
			copyState(AfterThetaR1, AfterThetaBack);
			copyState(BeforeChiR1, BeforeChiBack);
			copyState(BeforeThetaR2, TempBefTheta);
			copyState(AfterThetaR2, TempAftTheta);
			copyState(BeforeChiR2, TempBefChi);
			showCout();
		}
	}
	else {
		if (ComputeMinRevWeight(BeforeThetaBack) + ComputeWeight(BeforeChiBack) + ComputeWeight(TempBefChi) + ComputeWeight(TempLastBefChi) <= MaxCost4Round) {
			copyState(BeforeThetaR0, BeforeThetaBack);
			copyState(AfterThetaR0, AfterThetaBack);
			copyState(BeforeChiR0, BeforeChiBack);
			showCout4();

			copyState(BeforeThetaFinalR1, BeforeThetaBack);
			copyState(AfterThetaFinalR1, AfterThetaBack);
			copyState(BeforeChiFinalR1, BeforeChiBack);
			copyState(BeforeThetaFinalR2, TempBefTheta);
			copyState(AfterThetaFinalR2, TempAftTheta);
			copyState(BeforeChiFinalR2, TempBefChi);
			copyState(BeforeThetaFinalR3, BeforeTheta);
			copyState(AfterThetaFinalR3, AfterTheta);
			copyState(BeforeChiFinalR3, BeforeChi);
		}
	}
}

bool TrailSeedIterator::BackwardExtensionToChild() {
	if (C_Consistant()) {
		pos.push_back((pos[number - 1] - 1 + Size) % Size);
		number++;
		BeforeChiBack.flip(pos[number - 1]);
		BackwardUpdateState();
		return true;
	}
	else {
		return false;
	}
}

bool TrailSeedIterator::BackwardExtensionToSibling() {
	if (S_Consistant()) {
		BeforeChiBack.flip(pos[number - 1]);
		pos[number - 1] = (pos[number - 1] - 1 + Size) % Size;
		BeforeChiBack.flip(pos[number - 1]);
		BackwardUpdateState();
		return true;
	}
	else {
		return false;
	}
}

bool TrailSeedIterator::BackwardExtensionToPar() {
	if (number > 1) {
		BeforeChiBack.flip(pos[number - 1]);
		pos.pop_back();
		number--;
		return true;
	}
	else {
		return false;
	}
}

bool TrailSeedIterator::BackwardExtensionNext() {
	if (BackwardExtensionToChild() == false) {
		do {
			if (BackwardExtensionToSibling() == true) {
				return true;
			}
			else if (BackwardExtensionToPar() == false) {
				return false;
			}
		} while (true);
	}
	return true;
}

bool TrailSeedIterator::C_Consistant() {
	/* If the score exceeds the maximum weight.*/
	if (BackwardExtensionComputeLowerWeight(NumberOfRounds) == false) {
		return false;
	}
	/* If the weight exceeds the maximum weight or the position of the last bit reachs the end.*/
	if (ComputeWeight(BeforeChiBack) > maxBudgetBackward or pos[number - 1] == 0) {
		return false;
	}
	/* If input and output differences are not compatible.*/
	if (BeforeChiBack[Size - 2] == 0
		&& BeforeChiBack[Size - 1] == 0
		&& TempBefTheta[Size - 3] != BeforeChiBack[Size - 3]) {
		return false;
	}
	/* If the part of the input and output differences between two parameters are compatible.*/
	if (isCompatible(Size - 4, pos[number - 1])) {
		return true;
	}
	else {
		return false;
	}
}

bool TrailSeedIterator::S_Consistant() {
	/* If the score exceeds the maximum weight.*/
	if (BackwardExtensionComputeLowerWeight(NumberOfRounds) == false) {
		return false;
	}
	/* We say maxCost + 1, because sometimes going to sibling may decrease the weight,
	* for instance 10011 cost 4 but 10101 cost 3. */
	if (ComputeWeight(BeforeChiBack) > maxBudgetBackward + 1 or pos[number - 1] == 0) {
		return false;
	}
	/* If the part of the input and output differences between two parameters are compatible.*/
	if (isCompatible(Size - 4, pos[number - 1] + 1)) {
		return true;
	}
	else {
		return false;
	}
}

bool TrailSeedIterator::isCompatible(int aStart, int aEnd) {
	/* we need at least four bits of the input. */
	for (int i = aStart; i >= aEnd; i--) {
		if (BeforeChiBack[(i + 1) % Size] == 0
			&& BeforeChiBack[(i + 2) % Size] == 0
			&& TempBefTheta[i] != BeforeChiBack[i]) {
			return false;
		}
		else if (BeforeChiBack[i] == 1
			&& BeforeChiBack[(i + 1) % Size] == 1
			&& BeforeChiBack[(i + 2) % Size] == 0
			&& BeforeChiBack[(i + 3) % Size] == 1
			&& TempBefTheta[i] == TempBefTheta[(i + 1) % Size]) {
			return false;
		}
		else if (BeforeChiBack[i] == 0
			&& BeforeChiBack[(i + 1) % Size] == 1
			&& BeforeChiBack[(i + 2) % Size] == 0
			&& BeforeChiBack[(i + 3) % Size] == 1
			&& TempBefTheta[i] != TempBefTheta[(i + 1) % Size]) {
			return false;
		}
	}
	return true;
}

bool TrailSeedIterator::BackwardExtensionComputeLowerWeight(int aNumberOfRounds) {
	int ww = 0;
	for (int j = 0; j < number; j++) {
		int xx = 1;
		for (int i = 1; i < number; i++) {
			/*If the distance is less than (Size / piOffset)
			* then after passing the reverse of Pi map,
			* each piOffset bits at the state after Theta
			* has one active bit.
			* To cover this active bit, we need at least one active bit
			* at each 12 bits of the state before Theta.
			* Since each active bit before Theta covers the space of 8 bits,
			* then each active bit at state after Theta needs at least
			* one active bit at state before Theta that adds 2 to the weight. */
			if (pos[j] - pos[(j + i)% number] < Size / piOffset) {
				xx++;
				if (xx > ww) {
					ww = xx;
				}
			}
		}
	}
	if (aNumberOfRounds == 3) {
		if (Min_W0 + (2 * ww) + ComputeWeight(BeforeChiBack) + ComputeWeight(TempBefChi) > MaxCost4Round) {
			return false;
		}
	}
	else {
		if ((2 * ww) + ComputeWeight(BeforeChiBack) + ComputeWeight(TempBefChi) + ComputeWeight(TempLastBefChi) > MaxCost4Round) {
			return false;
		}
	}
	return true;
}
