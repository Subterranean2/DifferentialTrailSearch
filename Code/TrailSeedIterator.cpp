#include <iostream>
#include <stdio.h>
#include <fstream>
#include <iomanip>
#include<algorithm>
#include <bitset>
#include "TrailSeedIterator.h"
#include "Subterranean.h"

using namespace std;

TrailSeedIterator::TrailSeedIterator(int aForwardBackward, const int aMaxCost2Round, const int aMaxCost4Round
	, int his2Rout[], int his3Rout[], int his4Rout[]) :
	ForwardBackward(aForwardBackward), MaxCost2Round(aMaxCost2Round), MaxCost4Round(aMaxCost4Round) {
	for (int i = 0; i < 100; i++) {
		his2R[i] = his2Rout[i];
		his3R[i] = his3Rout[i];
		his4R[i] = his4Rout[i];
	}
	initialization();
	BeforeThetaR0.clear();
	AfterThetaR0.clear();
	BeforeChiR0.clear();
	BeforeThetaR1.clear();
	AfterThetaR1.clear();
	BeforeChiR1.clear();
	BeforeThetaR2.clear();
	AfterThetaR2.clear();
	BeforeChiR2.clear();
	BeforeThetaR3.clear();
	AfterThetaR3.clear();
	BeforeChiR3.clear();
	Three_RoundTrailCore(); 
}

void TrailSeedIterator::initialization() {
	GetParameters();
	InversePiParam();
	InverseThetaParam();
	for (int i = 0; i < Size; i++) {
		BeforeTheta.reset();
		AfterTheta.reset();
		BeforeChi.reset();
	}
	position.push_back(0);
	n = 1;
	updateState(position[n - 1]);

	if (ForwardBackward == 1) {
		makeHisto();
		if (computeRealWeight() <= MaxCost2Round) {
			his2R[computeRealWeight()]++;
		}
	}
}

bool TrailSeedIterator::toChild() {
	if (n > 1 && position[1] > Size / 2 + 1) {
		return false;
	}
	if (LowerBoundWeight() == true && position[n - 1] < Size - 1) {
		position.push_back(position[n - 1] + 1);
		n++;
		updateState(position[n - 1]);
		if (ForwardBackward == 1) {
			makeHisto();
			if (computeRealWeight() <= MaxCost2Round) {
				his2R[computeRealWeight()]++;
			}
		}
		return true;
	}
	else
		return false;
}

bool TrailSeedIterator::toSibling() {
	if (position[n - 1] < Size - 1 && LowerBoundWeight() == true && isCanonical() == true) {
		updateState(position[n - 1]);
		position[n - 1]++;
		updateState(position[n - 1]);
		if (ForwardBackward == 1) {
			makeHisto();
			if (computeRealWeight() <= MaxCost2Round) {
				his2R[computeRealWeight()]++;
			}
		}
		if (n == 3) {
			cout << "**New Node: ";
			for (int i = 0; i < Size; i++) {
				if (BeforeTheta[i] == 1)
					cout << i << ", ";
			}
			cout << "\n";
		}
		return true;
	}
	else
		return false;
}

bool TrailSeedIterator::toPar() {
	if (n > 2) {
		updateState(position[n - 1]); 
		position.pop_back();
		n--;
		return true;
	}
	else {
		return false;
	}
}


bool TrailSeedIterator::next() {
	if (toChild() == false) {
		do {
			if (toSibling() == true) {
				return true;
			}
			else if (toPar() == false) {
				return false;
			}
		} while (true);
	}
	else {
		return true;
	}
}

void TrailSeedIterator::updateState(int aP) {
	BeforeTheta.flip(aP);
	Theta(AfterTheta, aP);
	Pi(BeforeChi, aP);
}

int  TrailSeedIterator::ComputeLowerBoundMinRevWeight() {
	lowerBoundMinRevWeight = 0;
	/* To be able to compute a lower bound on the minimum reverse weight,
	* we need to compute the length of the non-zero string. */

	//Because the bit at position 0 before theta is always active the length is at least 1.
	int L = 1;
	for (int i = 0; i < Size; i++) { 
		if (BeforeTheta[i] == 1 && BeforeTheta[(i + 1) % Size] == 0
			&& BeforeTheta[(i + 2) % Size] == 0 && BeforeTheta[(i + 3) % Size] == 0
			&& BeforeTheta[(i + 4) % Size] == 0) {
			lowerBoundMinRevWeight += L / 6 + (L % 6 != 0) + L / 12;
			for (int k = i + 5; k < Size - 1; k++) {
				if (BeforeTheta[k] == 1) {
					L = 1;
					i = k - 1;
					k = Size;
				}
			}
		}
		else {
			L++;
		}
	}
	return lowerBoundMinRevWeight;
}

int  TrailSeedIterator::ComputeMinRevWeight(bitset<Size>& arr) {
	/* We make a bitset of a minimum reverse difference.*/
	bitset<Size> minRevDiff;

	int group1[] = { 141, 143, 145, 147, 157, 159,
		189, 191, 205, 207, 221, 223, 253, 255, 274, 278 };
	int inout;
	int index;
	bool fourZero = false;
	int startPoint;
	for (int i = Size - 1; i > 3; i--) {
		if (arr[i] == 0 && arr[i - 1] == 0 && arr[i - 2] == 0 && arr[i - 3] == 0) {
			startPoint = i - 4;
			fourZero = true;
			i = 0;
		}
	}
	if (fourZero == false) {
		return  Size / 6 + (Size % 6 != 0) + Size / 12 + (Size % 12 != 0);
	}

	/* Compute lower bound on the weight based on Section 5. */
	for (int i = Size + startPoint; i >= startPoint + 5; i--) {
		inout = 0;
		for (int y = Size + 2; y >= Size - 4; y--) {
			if (arr[(y + i) % Size] == 1) {
				inout += 1 << (Size + 2 - y);
			}
		}
		inout += 256 * minRevDiff[(i + 1) % Size] + 128 * minRevDiff[(i + 2) % Size];

		if (inout < 128) {//group1
			minRevDiff[(i) % Size] = arr[(i) % Size];
		}
		else {
			for (int j = 0; j < (sizeof(group1) / sizeof(group1[0])); j++) {//group1
				if (inout == group1[j]) {
					minRevDiff[(i) % Size] = 1;
				}
			}

			if (inout == 209 or inout == 211 or inout == 342 or inout == 338) {//group2
				bool odd = 0;
				for (int xx = 5; xx < Size; xx += 2) {
					index = i - xx + Size;
					if (arr[index % Size] == 0 && arr[(index - 1) % Size] == 1) {
						odd ^= 1;
						if (arr[(index - 2) % Size] == 0 && arr[(index - 3) % Size] == 0 && odd) {
							minRevDiff[(i) % Size] = 1;
						}
					}
					else {
						break;
					}
				}
			}
			if (inout == 173 or inout == 175) {//group2
				bool odd = 0;
				for (int xx = 5; xx < Size; xx += 2) {
					index = i - xx + Size;
					if (arr[index % Size] == 1 && arr[(index - 1) % Size] == 0) {
						odd ^= 1;
						if (arr[(index - 2) % Size] == 0 && odd) {
							minRevDiff[(i) % Size] = 1;
						}
					}
					else {
						break;
					}
				}
			}
		}
	}
	return ComputeWeight(minRevDiff);
}

int  TrailSeedIterator::ComputeLowerBoundBefChiWeight() {
	/* We make a bitset of all active bits that count.
	This weight counts and will not goes away because there will be no overlap. */
	bitset<Size> BeforeChilowerBoundBefChiWeight;
	/* If the position of an active bit is less than (position[n - 1] - thetaOffset.back())
	then it counts. */
	for (unsigned int i = Size - thetaOffset.back(); i < position[n - 1] + Size - thetaOffset.back(); i++) {
		if (AfterTheta[i % Size] == 1) {
			BeforeChilowerBoundBefChiWeight.flip((i * piOffsetInverse) % Size);
		}
	}
	if (position[n - 1] >= Size - thetaOffset.back()) {
		cout << "There is an active bit in the last 8 bits.";
	}
	int addedWeightByLastActiveBit = 0;
	if (n > 1) {
		if (position[n - 1] - position[n - 2] > (2 * piOffset + thetaOffset.back()) && position[n - 1] - position[0] < Size - (2 * piOffset + thetaOffset.back())) {
			addedWeightByLastActiveBit = 6;
		}
		else if (position[n - 1] - position[n - 2] > (2 * piOffset + thetaOffset[1]) && position[n - 1] - position[0] < Size - (2 * piOffset + thetaOffset[1])) {
			addedWeightByLastActiveBit = 5;
		}
		else if (position[n - 1] - position[n - 2] > (2 * piOffset + thetaOffset[0]) && position[n - 1] - position[0] < Size - (2 * piOffset + thetaOffset[0])) {
			addedWeightByLastActiveBit = 4;
		}
		else if (position[n - 1] - position[n - 2] > thetaOffset.back() && position[n - 1] - position[0] < thetaOffset.back()) {
			addedWeightByLastActiveBit = 3;
		}
	}
	/* the weight of b' is always greater than 5 in any 3-round trail core with weight smaller than 61. */
	return max(ComputeWeight(BeforeChilowerBoundBefChiWeight) + addedWeightByLastActiveBit, 5);
}

int TrailSeedIterator::ComputeWeight(bitset<Size> &arr) {
	int Weightt = 0;
	for (int i = 0; i < Size; i++) {
		if (arr[i] == 1) {
			Weightt++;
			if (arr[(i + 1) % Size] == 0 && arr[(i + 2) % Size] == 0) {
				Weightt++;
			}
		}
	}
	return Weightt;
}

bool TrailSeedIterator::LowerBoundWeight() {
	if (n > 4 && position[n - 4] == Size - 4) {
		return false;
		cout << "Fully active state before theta!\n";
	}
	if (ComputeLowerBoundMinRevWeight() + ComputeLowerBoundBefChiWeight() > MaxCost2Round) {
		return false;
	}
	return true;
}

int TrailSeedIterator::computeRealWeight() {
	if (n > 4 && position[n - 4] == Size - 4){
		return false;
		cout << "Fully active state before theta!\n";
	}
	return (ComputeMinRevWeight(BeforeTheta) + ComputeWeight(BeforeChi));
}

void TrailSeedIterator::makeHisto() {
	if (LowerBoundWeight() == true && computeRealWeight() <= MaxCost2Round) {
		his[ComputeMinRevWeight(BeforeTheta)][ComputeWeight(BeforeChi)]++;
		hisss[ComputeMinRevWeight(BeforeTheta) + ComputeWeight(BeforeChi)]++;
	}
}

bool TrailSeedIterator::isCanonical() {
	if (n == 1 && position[n - 1] != 0) {
		return false;
	}
	for (int i = n - 1; i > 0; i--) {
		int h = 1;
		unsigned int x = Size - position[n - 1];
		for (int j = n + i; j > i + 1; j--) {
			if (x < (Size + position[j % n] - position[(j - 1) % n]) % Size) {
				return false;
			}
			else if (x == (Size + position[j % n] - position[(j - 1) % n]) % Size) {
				h++;
				if (h == n) {
					 cout << "Symmetry!!!!!!!!!\nThe only possible symmetry is when we have fully active state!";
				}
				else {
					x = position[n - h + 1] - position[n - h];
				}
			}
			else {
				j = i;
			}
		}
	}
	return true;
}

void TrailSeedIterator::hi22(ofstream& fou) {
	fou << "   2-round \n";
	fou << setw(4) << "w(a)/w(b)" << ",";
	for (int i = 2; i <= MaxCost2Round - 2; i++) {
		fou << setw(3) << i << ",";
	}
	fou << "\n";
	for (int i = 2; i <= MaxCost2Round - 2; i++) {
		fou << setw(4) << i << ",     ";
		for (int j = 2; j <= MaxCost2Round - i; j++) {
			fou << setw(3) << his[i][j] << ",";
		}
		fou << "\n";
	}
	fou << "\n\n\n\n\n\n\n";
	for (int i = 8; i <= MaxCost2Round; i++) {
		fou << setw(4) << i << ", " << hisss[i] << "\n";
	}
}

void TrailSeedIterator::hi2(ofstream& fou) {
	fou << "   2-round \n";
	fou << setw(4) << "w" << setw(16) << "#";
	fou << "\n";
	for (int i = 8; i <= MaxCost2Round; i++) {
		fou << setw(4) << i;
		if (his2R[i] == 0)
			fou << setw(16) << "-";
		else
			fou << setw(16) << his2R[i];
		fou << "\n";
	}
}

void TrailSeedIterator::hi(ofstream& fou) {
	fou << "   3-round \n";
	fou << setw(4) << "w" << setw(16) << "#";
	fou << "\n";
	for (int i = 25; i <= MaxCost4Round - Min_W2; i++) {
		fou << setw(4) << i;
		if (his3R[i] == 0)
			fou << setw(16) << "-";
		else
			fou << setw(16) << his3R[i];
		fou << "\n";
	}
}

void TrailSeedIterator::disssplay(vector<bool>& arr, int j, ofstream& fou) {
	for (int i = 0; i < Size; i++) {
		if (arr[j * Size + i] == 0)
			fou << ".";
		if (arr[j * Size + i] == 1)
			fou << "1";
	}
}

void TrailSeedIterator::displayCout(vector<bool>& arr, int j) {
	for (int i = 0; i < Size; i++) {
		if (arr[j * Size + i] == 0)
			cout << ".";
		if (arr[j * Size + i] == 1)
			cout << "1";
	}
}

void TrailSeedIterator::showCout() {
	if (BeforeThetaR1.size() > 0) {
		int j = (BeforeThetaR1.size() / Size) - 1;
		cout << "\n\n\n\nw_rev(a'0) + w(b'0) + w(b'1) = ";
		if (ForwardBackward == 0) {//backward
			cout << ComputeMinRevWeight(BeforeThetaBack);
			cout << " + " << ComputeWeight(BeforeChiBack);
			cout << " + " << ComputeWeight(BeforeChi);
			cout << " = " << ComputeMinRevWeight(BeforeThetaBack) + ComputeWeight(BeforeChiBack) + ComputeWeight(BeforeChi);
		}
		if (ForwardBackward == 1) {//forward
			bitset<Size> xRb1;
			if (BeforeThetaR2.size() > 1) {
				for (int i = 0; i < Size; i++) {
					xRb1[i] = BeforeChiR2[BeforeThetaR2.size() - Size + i];
				}
			}
			cout << ComputeMinRevWeight(BeforeTheta);
			cout << " + " << ComputeWeight(BeforeChi);
			cout << " + " << ComputeWeight(xRb1);
			cout << " = " << ComputeMinRevWeight(BeforeTheta) + ComputeWeight(BeforeChi) + ComputeWeight(xRb1);
		}
		cout << "\na1'    " << setw(3);
		displayCout(BeforeThetaR1, j);
		cout << "\nAf.Th. " << setw(3);
		displayCout(AfterThetaR1, j);
		cout << "\nb1'    " << setw(3);
		displayCout(BeforeChiR1, j);
		cout << "\na2'    " << setw(3);
		displayCout(BeforeThetaR2, j);
		cout << "\nAf.Th. " << setw(3);
		displayCout(AfterThetaR2, j);
		cout << "\nb2'    " << setw(3);
		displayCout(BeforeChiR2, j);
	}
}

void TrailSeedIterator::showCout4() {
	if (ForwardBackward == 0) {
		cout << "\n\n\n\n\n\n Number of four-round trails BACKEARD = " << (BeforeThetaR0.size() / Size);
		if (BeforeThetaR0.size() > 0) {
			int j = (BeforeThetaR0.size() / Size) - 1;
			cout << "\n\n\nw_rev(a'0) + w(b'0) + w(b'1) + w(b'2) = ";
			cout << ComputeMinRevWeight(BeforeThetaBack);
			cout << " + " << ComputeWeight(BeforeChiBack);
			cout << " + " << ComputeWeight(TempBefChi);
			cout << " + " << ComputeWeight(TempLastBefChi);
			cout << " = " << ComputeMinRevWeight(BeforeThetaBack) + ComputeWeight(BeforeChiBack)
				+ ComputeWeight(TempBefChi) + ComputeWeight(TempLastBefChi);

			cout << "\na1'    " << setw(3);
			displayCout(BeforeThetaR0, j);
			cout << "\nAf.Th. " << setw(3);
			displayCout(AfterThetaR0, j);
			cout << "\nb1'    " << setw(3);
			displayCout(BeforeChiR0, j);
			cout << "\na2'    " << setw(3);
			displayCout(BeforeThetaR1, j);
			cout << "\nAf.Th. " << setw(3);
			displayCout(AfterThetaR1, j);
			cout << "\nb2'    " << setw(3);
			displayCout(BeforeChiR1, j);
			cout << "\na3'    " << setw(3);
			displayCout(BeforeThetaR2, j);
			cout << "\nAf.Th. " << setw(3);
			displayCout(AfterThetaR2, j);
			cout << "\nb3'    " << setw(3);
			displayCout(BeforeChiR2, j);
			cout << "\n\n\n\n";
		}
	}

	if (ForwardBackward == 1) {
		cout << "\n\n\n\n\n\n Number of four-round trails FORWARD = " << (BeforeThetaR3.size() / Size) << "\n\n";
		if (BeforeThetaR3.size() > 0) {
			int j = (BeforeThetaR3.size() / Size) - 1;
			cout << "\n\n\nw_rev(a'0) + w(b'0) + w(b'1) + w(b'2) = ";
			cout << ComputeMinRevWeight(BeforeTheta);
			cout << " + " << ComputeWeight(BeforeChi);
			cout << " + " << ComputeWeight(TempBefChi);
			cout << " + " << ComputeWeight(OffsetBeforeChiF);
			cout << " = " << ComputeMinRevWeight(BeforeTheta) + ComputeWeight(BeforeChi)
				+ ComputeWeight(TempBefChi) + ComputeWeight(OffsetBeforeChiF);
			cout << "\na1'    " << setw(3);
			displayCout(BeforeThetaR1, j);
			cout << "\nAf.Th. " << setw(3);
			displayCout(AfterThetaR1, j);
			cout << "\nb1'    " << setw(3);
			displayCout(BeforeChiR1, j);
			cout << "\na2'    " << setw(3);
			displayCout(BeforeThetaR2, j);
			cout << "\nAf.Th. " << setw(3);
			displayCout(AfterThetaR2, j);
			cout << "\nb2'    " << setw(3);
			displayCout(BeforeChiR2, j);
			cout << "\na3'    " << setw(3);
			displayCout(BeforeThetaR3, j);
			cout << "\nAf.Th. " << setw(3);
			displayCout(AfterThetaR3, j);
			cout << "\nb3'    " << setw(3);
			displayCout(BeforeChiR3, j);
			cout << "\n\n\n\n";
		}
	}
}

void TrailSeedIterator::showw4(ofstream& fout) {
	if (ForwardBackward == 0) {
		fout << "Number of four-round trail cores BACKWARD = " << (BeforeThetaFinalR1.size() / Size) << "\n\n";
		if (BeforeThetaFinalR1.size() > 0) {
			bitset<Size> xRa1;
			bitset<Size> xRb1;
			bitset<Size> xRb2;
			bitset<Size> xRb3;
			for (unsigned int j = 0; j < BeforeThetaFinalR1.size() / Size; j++) {
				for (int i = 0; i < Size; i++) {
					xRa1[i] = BeforeThetaFinalR1[j * Size + i];
					xRb1[i] = BeforeChiFinalR1[j * Size + i];
					xRb2[i] = BeforeChiFinalR2[j * Size + i];
					xRb3[i] = BeforeChiFinalR3[j * Size + i];
				}
				fout << "\n\n\nw_rev(a1') + w(b1') + w(b2') + w(b3') = ";
				fout << ComputeMinRevWeight(xRa1);
				fout << " + " << ComputeWeight(xRb1);
				fout << " + " << ComputeWeight(xRb2);
				fout << " + " << ComputeWeight(xRb3);
				fout << " = " << ComputeMinRevWeight(xRa1) + ComputeWeight(xRb1)
					+ ComputeWeight(xRb2) + ComputeWeight(xRb3);
				his4R[ComputeMinRevWeight(xRa1) + ComputeWeight(xRb1) + ComputeWeight(xRb2) + ComputeWeight(xRb3)]++;
				fout << "\na1'    " << setw(3);
				disssplay(BeforeThetaFinalR1, j, fout);
				fout << "\nAf.Th. " << setw(3);
				disssplay(AfterThetaFinalR1, j, fout);
				fout << "\nb1'    " << setw(3);
				disssplay(BeforeChiFinalR1, j, fout);
				fout << "\na2'    " << setw(3);
				disssplay(BeforeThetaFinalR2, j, fout);
				fout << "\nAf.Th. " << setw(3);
				disssplay(AfterThetaFinalR2, j, fout);
				fout << "\nb2'    " << setw(3);
				disssplay(BeforeChiFinalR2, j, fout);
				fout << "\na3'    " << setw(3);
				disssplay(BeforeThetaFinalR3, j, fout);
				fout << "\nAf.Th. " << setw(3);
				disssplay(AfterThetaFinalR3, j, fout);
				fout << "\nb3'    " << setw(3);
				disssplay(BeforeChiFinalR3, j, fout);
				fout << "\n\n\n\n";
			}
		}
	}

	if (ForwardBackward == 1) {
		fout << "\n\n\n\n\n\n\n\n\n\n\n";
		fout << "_________________________________________________________________________________________\n";
		fout << "Number of four - round trail cores FORWARD = " << (BeforeThetaFinalR3.size() / Size) << "\n\n";
		if (BeforeThetaFinalR3.size() > 0) {
			bitset<Size> xRa1;
			bitset<Size> xRb1;
			bitset<Size> xRb2;
			bitset<Size> xRb3;
			for (unsigned int j = 0; j < BeforeThetaFinalR3.size() / Size; j++) {
				for (int i = 0; i < Size; i++) {
					xRa1[i] = BeforeThetaFinalR1[j * Size + i];
					xRb1[i] = BeforeChiFinalR1[j * Size + i];
					xRb2[i] = BeforeChiFinalR2[j * Size + i];
					xRb3[i] = BeforeChiFinalR3[j * Size + i];
				}
				fout << "\n\n\nw_rev(a1') + w(b1') + w(b2') + w(b3') = ";
				fout << ComputeMinRevWeight(xRa1);
				fout << " + " << ComputeWeight(xRb1);
				fout << " + " << ComputeWeight(xRb2);
				fout << " + " << ComputeWeight(xRb3);
				fout << " = " << ComputeMinRevWeight(xRa1) + ComputeWeight(xRb1)
					+ ComputeWeight(xRb2) + ComputeWeight(xRb3);
				his4R[ComputeMinRevWeight(xRa1) + ComputeWeight(xRb1) + ComputeWeight(xRb2) + ComputeWeight(xRb3)]++;
				fout << "\na1'    " << setw(3);
				disssplay(BeforeThetaFinalR1, j, fout);
				fout << "\nAf.Th. " << setw(3);
				disssplay(AfterThetaFinalR1, j, fout);
				fout << "\nb1'    " << setw(3);
				disssplay(BeforeChiFinalR1, j, fout);
				fout << "\na2'    " << setw(3);
				disssplay(BeforeThetaFinalR2, j, fout);
				fout << "\nAf.Th. " << setw(3);
				disssplay(AfterThetaFinalR2, j, fout);
				fout << "\nb2'    " << setw(3);
				disssplay(BeforeChiFinalR2, j, fout);
				fout << "\na3'    " << setw(3);
				disssplay(BeforeThetaFinalR3, j, fout);
				fout << "\nAf.Th. " << setw(3);
				disssplay(AfterThetaFinalR3, j, fout);
				fout << "\nb3'    " << setw(3);
				disssplay(BeforeChiFinalR3, j, fout);
				fout << "\n\n\n\n\n\n\n\n\n";
			}
		}
		fout << "\n_________________________________________________________________________________________\n";
	}
}

void TrailSeedIterator::showw3(ofstream& fout) {
	if (ForwardBackward == 0) {
		fout << "Number of three-round trail cores BACKWARD = ";
		fout << (BeforeThetaR1.size() / Size) << "\n\n";
	}
	if (ForwardBackward == 1) {
		fout << "\n\n\n\n\n\n\n\n\n\n\n\n";
		fout << "_________________________________________________________________________________________\n";
		fout << "Number of three-round trail cores FORWARD = ";
		fout << (BeforeChiR2.size() / Size) << "\n\n";
	}
	bitset<Size> xRa1;
	bitset<Size> xRb1;
	bitset<Size> xRb2;
	for (unsigned int j = 0; j < BeforeThetaR1.size() / Size; j++) {
		for (int i = 0; i < Size; i++) {
			xRa1[i] = BeforeThetaR1[j * Size + i];
			xRb1[i] = BeforeChiR1[j * Size + i];
			xRb2[i] = BeforeChiR2[j * Size + i];
		}
		fout << "\n\n\n\nw_rev(a'0) + w(b'0) + w(b'1) = ";
		fout << ComputeMinRevWeight(xRa1);
		fout << " + " << ComputeWeight(xRb1);
		fout << " + " << ComputeWeight(xRb2);
		fout << " = " << ComputeMinRevWeight(xRa1) + ComputeWeight(xRb1) + ComputeWeight(xRb2);
		fout << "\na'0    " << setw(3);
		disssplay(BeforeThetaR1, j, fout);
		fout << "\nAf.Th. " << setw(3);
		disssplay(AfterThetaR1, j, fout);
		fout << "\nb'0    " << setw(3);
		disssplay(BeforeChiR1, j, fout);
		fout << "\na'1    " << setw(3);
		disssplay(BeforeThetaR2, j, fout);
		fout << "\nAf.Th. " << setw(3);
		disssplay(AfterThetaR2, j, fout);
		fout << "\nb'1    " << setw(3);
		disssplay(BeforeChiR2, j, fout);
		if (ForwardBackward == 0) {
			his3R[ComputeMinRevWeight(xRa1) + ComputeWeight(xRb1) + ComputeWeight(xRb2)]++;
		}
		if (ForwardBackward == 1 && ComputeWeight(xRb1) + ComputeWeight(xRb2) >= MaxCost4Round - MaxCost2Round ) {
			his3R[ComputeMinRevWeight(xRa1) + ComputeWeight(xRb1) + ComputeWeight(xRb2)]++;
		}
	}
}


void TrailSeedIterator::machinReadable(ofstream& fout) {
	int s1, s2;
	for (unsigned int j = 0; j < BeforeThetaR1.size() / Size; j++) {
		for (int i = 0; i < Size; i++) {
			if (BeforeThetaR1[j * Size + i] == 1) {
				fout << i;
				s1 = i+1;
				i = Size;
			}
		}
		for (int i = s1; i < Size; i++) {
			if (BeforeThetaR1[j * Size + i] == 1) {
				fout << ", " << i;
			}
		}
		fout << "\n";
		for (int i = 0; i < Size; i++) {
			if (BeforeThetaR2[j * Size + i] == 1) {
				fout << i;
				s2 = i+1;
				i = Size;
			}
		}
		for (int i = s2; i < Size; i++) {
			if (BeforeThetaR2[j * Size + i] == 1) {
				fout << ", " << i;
			}
		}
		fout << "\n\n";
	}
}
