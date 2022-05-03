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

void TrailSeedIterator::Three_RoundTrailCore() {

	if (isCanonical() == true) {
		/* Backward extension */
		if (ForwardBackward == 0 && ComputeMinRevWeight(BeforeTheta) + ComputeWeight(BeforeChi) <= MaxCost2Round) {
			if (ComputeMinRevWeight(BeforeTheta) + ComputeWeight(BeforeChi) <= MaxCost4Round - MinWeight2_Round) {//check me
				/* We need to generate all 4-r trail cores up to MaxCost4Round.
				* It means w_rev(a'0) + w(b'1) + w(b'2) + w(b'3) <= MaxCost4Round.
				* min(w_rev(a'0)+ w(b'1)) = 8, and w(b'3) is known.
				* Therefore, the maximum of w(b'1) can be computed as follows: */
				maxBudgetBackward = MaxCost4Round - MinWeight2_Round - ComputeWeight(BeforeChi);
				/* We know w(b'2) + w(b'3) <= MaxCost2Round
				* so, (w(b'3) <= MaxCost2Round - w(b'3) < maxBudgetBackward)*/
				if (MaxCost2Round - ComputeWeight(BeforeChi) < maxBudgetBackward) {
					maxBudgetBackward = MaxCost2Round - ComputeWeight(BeforeChi);
				}
				/* We know w(b'1) >= 5. */
				for (int j = 0; j < Size; j++) {
					TempBefTheta[j] = BeforeTheta[j];
					TempAftTheta[j] = AfterTheta[j];
					TempBefChi[j] = BeforeChi[j];
				}
				NumberOfRounds = 3;
				BackwardExtension();
			}
		}

		/* Forward extension */
		if (ForwardBackward == 1 && ComputeMinRevWeight(BeforeTheta) + ComputeWeight(BeforeChi) <= MaxCost2Round) {
			/*We need to generate all 4-r trail cores up to MaxCost4Round.
			* It means w_rev(a'0) + w(b'1) + w(b'2) + w(b'3) <= MaxCost4Round.
			* we know min(w(b'3)) = 5 therefore: */
			if (ComputeMinRevWeight(BeforeTheta) + ComputeWeight(BeforeChi) <= MaxCost4Round - Min_W2) {
				for (int j = 0; j < Size; j++) {
					TempBefTheta[j] = BeforeTheta[j];
					TempAftTheta[j] = AfterTheta[j];
					TempBefChi[j] = BeforeChi[j];
				}
				NumberOfRounds = 3;
				ForwardExtension();
			}

		}
	}
}

void TrailSeedIterator::copyState(vector<bool>& arr1, bitset<Size>& arr2) {
	for (int i = 0; i < Size; i++) {
		arr1.push_back(arr2[i]);
	}
}	