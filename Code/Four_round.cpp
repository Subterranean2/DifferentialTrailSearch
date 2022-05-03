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

void TrailSeedIterator::Four_RoundTrailCore() {
	/* Backward extension */
	if (ForwardBackward == 0) {
		for (unsigned int i = 0; i < BeforeThetaR1.size() / Size; i++) {
			for (int j = 0; j < Size; j++) {
				TempBefTheta[j] = BeforeThetaR1[i * Size + j];
				TempAftTheta[j] = AfterThetaR1[i * Size + j];
				TempBefChi[j] = BeforeChiR1[i * Size + j];

				TempLastBefChi[j] = BeforeChiR2[i * Size + j];
			}
			/*We need to generate all 4-r trail cores up to MaxCost4Round.
			* It means w_rev(a'0) + w(b'1) + w(b'2) + w(b'3) <= MaxCost4Round,
			* min(w_rev(a'0)) = 2 and w(b'2) and w(b'3) are known.
			* Therefore, maximum of w(b'1) can be computed as follows:*/
			maxBudgetBackward = MaxCost4Round - Min_W0 - ComputeWeight(TempBefChi) - ComputeWeight(TempLastBefChi);
			if (ComputeMinRevWeight(TempBefTheta) <= maxBudgetBackward) {
				NumberOfRounds = 4;
				BackwardExtension();
				cout << "   " << i;
			}
		}
	}

	/* Forward extension */
	if (ForwardBackward == 1) {
		for (unsigned int i = 0; i < BeforeThetaR2.size() / Size; i++) {
			for (int j = 0; j < Size; j++) {
				TempBefTheta[j] = BeforeThetaR2[i * Size + j];
				TempAftTheta[j] = AfterThetaR2[i * Size + j];
				TempBefChi[j] = BeforeChiR2[i * Size + j];
			}
			NumberOfRounds = 4;
			ForwardExtension();
			cout << "   " << i;
		}
		cout << "\n\n\n\n\n\n #(four-round trails Forward extension) = " << BeforeThetaR3.size() / Size;
	}
}