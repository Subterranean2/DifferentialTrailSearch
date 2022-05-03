#pragma once
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <algorithm>
#include <functional>
#include <array>
#include <vector> 
#include "Subterranean.h"


using namespace std;

/* Class TrailSeedIterator represent an iterator for valid (orbital) point in tree*/
class TrailSeedIterator : public Subterranean
{
public:
	/* Number of active points. */
	int n;

	/* Position of active points. */
	vector<unsigned int> position;

	/* This parameter represent the maximum weight of 2 rounds. */
	int MaxCost2Round;

	/* This parameter represent the maximum weight of 3 rounds. */
	int MaxCost3Round;

	/* This parameter represent the maximum weight of 4 rounds. */
	int MaxCost4Round;

	int MinWeightAfterChi;
	int WeightBeforeChi;
public:
	/* TrailSeedIterator constructor.
	* The first parameter represents the direction of extension,
	* the second parameter represents the maximum weight of 2-round trail cores, 
	* the third parameter represents the maximum weight of 4-round trail cores, 
	* the remaining parameters represent the arrays for storing the histograms. */
	TrailSeedIterator(int aForwardBackward, const int aMaxCost2Round, const int aMaxCost4Round
		, int his2Rout[], int his3Rout[], int his4Rout[]);

	/* Method for initialization.
	* It adds one active bit to the root of the tree. */
	void initialization();

	/* Method to update the whole trial after adding each active bit. */
	void updateState(int aP);

	/* Method to go to the next point, or returns false if no point is acceptable. */
	bool next();

	/* Method to compute the score of a'. */
	int ComputeLowerBoundMinRevWeight();

	/* Parameter to store the score of a'. */
	int lowerBoundMinRevWeight;

	/* Method to compute the minimum reverse weight of a'. */
	int ComputeMinRevWeight(bitset<Size>& arr);

	/* Method to compute the score of b'. */
	int ComputeLowerBoundBefChiWeight();

	/* Method to the weight of a given state before Chi. */
	int ComputeWeight(bitset<Size>& arr);

	/* Method to compute the weight of a given 2-round trail core. */
	int computeRealWeight();

	bool LowerBoundWeight();

	/* Array to store the histogram of all 2-round trail cores. */
	int his2R[100] = {};

	/* Array to store the histogram of all 3-round trail cores. */
	int his3R[100] = {};

	/* Array to store the histogram of all 4-round trail cores. */
	int his4R[100] = {};

	/* Array to store the histogram of all 2-round trail cores
	* regarding the values of w_rev(a'0) and w(b'1). */
	int his[70][70] = { {} };

	/* Array to store the histogram of all 2-round trail cores. */
	int hisss[50] = {};

	/* Method to generate a histogram of all 2-round trail cores. */
	void makeHisto();

	/* Method to represent number of 3-round trail cores. */
	void hi(ofstream& fou);

	/* Method to represent number of 2-round trail cores. */
	void hi2(ofstream& fou);

	/* Method to represent number of 2-round trail cores in details. */
	void hi22(ofstream& fou);

	/* Method to represent 3-round trail cores. */
	void showw3(ofstream& fou);

	/* Method to represent 4-round trail cores. */
	void showw4(ofstream& fou);

	/* Method to represent an intermediate state. */
	void disssplay(vector<bool>& arr, int j, ofstream& fou);

	/* Method to generate a machine readable file of required 3-round trail cores. */
	void machinReadable(ofstream& fou);

	/* Method to represent 3-round trail cores in the output console. */
	void showCout();

	/* Method to represent 4-round trail cores in the output console. */
	void showCout4();

	/* Method to represent an intermediate state in the output console. */
	void displayCout(vector<bool>& arr, int j);

protected:
	/* Method to go to the first child.
	* It returns false if we exceed maxCost or reach last available point.
	* Otherwise, returns true and add a new point to tree. */
	bool toChild();

	/* Method to go to the sibling.
	* It doesn't add a new point to tree and just consider other sibling of current point.
	* It returns false if there is no sibling available.
	* Otherwise, returns true and change last point of unit-list. */
	bool toSibling();

	/* Method to go back to the parent node.
	* This method discard last point of unit-list. */
	bool toPar();

	/** Parameter to show if a state is canonical or not. */
	bool isCanonical();

public:
	/* Method to generate required 3-round trail cores. */
	void Three_RoundTrailCore();

	/* Parameter to specify the direction of extension.
	* 0 represents backward extension and 1 represents forward. */
	int ForwardBackward = -1;

	/* Method to generate all 4-round trail cores. */
	void Four_RoundTrailCore();

	/* Method to extend trail cores in the forward direction. */
	void ForwardExtension();

	/* Method to generate bases vectors with the minimum Hamming weight. */
	void GenerateBases(int aType, int aP);

	/* Method to form affine space by adding basis vectors to the offset. */
	void updateOffset(int aP);

	/* Method to generate all basis vectors with minimum Hamming weihgt and canonical offset. */
	void setBasisOffset();

	/* Parameter to store the state before Theta. */
	bitset<Size> OffsetBeforeThetaF;

	/* Parameter to store the Basis vectors before Theta. */
	vector<unsigned int> BasisBeforeThetaF;

	/* Parameter to store the state after Theta. */
	bitset<Size> OffsetAfterThetaF;

	/* Parameter to store the Basis vectors after Theta. */
	vector<unsigned int> BasisAfterThetaF;

	/* Parameter to store the state before Chi. */
	bitset<Size> OffsetBeforeChiF;

	/* Parameter to store the Basis vectors before Chi. */
	vector<unsigned int> BasisBeforeChiF;

	/* Parameter to store state of the union of all possible active bits of offset and bases. */
	vector<unsigned int> allOutputActiveBits;

	/* Parameter that represents the minimum weight of the last round
	* in any 3-round trial core with weight smaller than 60. */
	int Min_W2 = 5;

	/* Parameter that represents the lowest minimum reverse weight of the first round
	* in any 3-round trial core with weight smaller than 60. */
	int Min_W0 = 4;

	/* Parameter that represents the minimum weight of any 2-round trial core. */
	int MinWeight2_Round = 8;

public:
	/* Method to extend trails in the backward direction. */
	void BackwardExtension();
	/* Method to generate the first child of a node in the tree.
	* It returns false if it is not "C_Consistant()",
	* otherwise, returns true and add a new node to tree. */
	bool BackwardExtensionToChild();

	/* Method to generate the first sibling of a node in the tree.
	* It returns false if it is not "S_Consistant()",
	* otherwise, returns true and add a new node to tree. */
	bool BackwardExtensionToSibling();

	/* Method to go back to the parent node.
	* It returns false if no active bit remains,
	* otherwise, returns true and pops the last unit form the unit-list. */
	bool BackwardExtensionToPar();

	/* Method to go to the next input difference for a given output difference, 
	* or returns false if no point is acceptable. */
	bool BackwardExtensionNext();

	/* Method to initialize backward extension. */
	void GenerateAllInputDiffInitialize();

	/* Method to check if a node is eligible to have a child. */
	bool C_Consistant();

	/* Method to check if a node is eligible to have a sibling. */
	bool S_Consistant();

	/* Method to compute the lower bound on the weight during backward extension. */
	bool BackwardExtensionComputeLowerWeight(int aNumberOfRounds);

	/* Method to compute the last round's score. */
	bool ComputeLastRoundScore(int aNumberOfRounds);

	/* Method to update the whole trail during backward extension. */
	void BackwardUpdateState();

	/* Parameter to store the score of the extended round. */
	int ww = 0;

	/* Parameter to store the number of active bits of the input difference
	* for a given output difference. */
	int number;

	/* Parameter to store the maximum weight of the input difference
	* for a given output difference. */
	int maxBudgetBackward;

	/* Parameter to store the position of active bits of the input difference
	* for a given output difference. */
	vector<unsigned int> pos;

	/* Parameter to store the value of the state before Theta
	* in the backward extended round. */
	bitset<Size> BeforeThetaBack;

	/* Parameter to store the value of the state after Theta
	* in the backward extended round. */
	bitset<Size> AfterThetaBack;

	/* Parameter to store the value of the state before Chi
	* in the backward extended round. */
	bitset<Size> BeforeChiBack;

	/* Parameter to temporarily store the value of the state before Theta. */
	bitset<Size> TempBefTheta;

	/* Parameter to temporarily store the value of the state after Theta. */
	bitset<Size> TempAftTheta;

	/* Parameter to temporarily store the value of the state before Chi. */
	bitset<Size> TempBefChi;

	/* Parameter to temporarily store the value of the state before Chi
	* in the last round. */
	bitset<Size> TempLastBefChi;

public:
	/* Method to copy a state. */
	void copyState(vector<bool>& arr1, bitset<Size>& arr2);

	/* Method to copy a state in the forward extension. */
	void copyStateForward(int aNumberOfRounds);

	/* Method to copy a state in the backward extension. */
	void copyStateBackward(int aNumberOfRounds);

	/* Parameter to store the number of rounds. */
	int NumberOfRounds;

	/* Method to check if input and output differences are partialy compatible
	* between positions "aStart" and "aEnd". */
	bool isCompatible(int aStart, int aEnd);

	/* Parameter to represent the value of the state
	* before Theta of the first round. */
	vector<bool> BeforeThetaFinalR1;

	/* Parameter to represent the value of the state
	* after Theta of the first round. */
	vector<bool> AfterThetaFinalR1;

	/* Parameter to represent the value of the state
	* before Chi of the first round. */
	vector<bool> BeforeChiFinalR1;

	/* Parameter to represent the value of the state
	* before Theta of the second round. */
	vector<bool> BeforeThetaFinalR2;

	/* Parameter to represent the value of the state
	* after Theta of the second round. */
	vector<bool> AfterThetaFinalR2;

	/* Parameter to represent the value of the state
	* before Chi of the second round. */
	vector<bool> BeforeChiFinalR2;

	/* Parameter to represent the value of the state
	* before Theta of the third round. */
	vector<bool> BeforeThetaFinalR3;

	/* Parameter to represent the value of the state
	* after Theta of the third round. */
	vector<bool> AfterThetaFinalR3;

	/* Parameter to represent the value of the state
	* before Chi of the third round. */
	vector<bool> BeforeChiFinalR3;
public:
	/* Parameter BeforeThetaR0 represents the value of the state
	* before Theta of round 0 (after Chi). */
	vector<bool> BeforeThetaR0;

	/* Vector AfterThetaR0 represents the value of the state
	* after Theta of round 0. */
	vector<bool> AfterThetaR0;

	/* Parameter BeforeChiR0 represents the value of the state
	* before Chi of round 0 (after Pi). */
	vector<bool> BeforeChiR0;

	/* Parameter BeforeThetaR1 represents the value of the state
	* before Theta of round 1 (after Chi). */
	vector<bool> BeforeThetaR1;

	/* Vector AfterThetaR1 represents the value of the state
	* after Theta of round 1. */
	vector<bool> AfterThetaR1;

	/* Parameter BeforeChiR1 represents the value of the state
	* before Chi of round 1 (after Pi). */
	vector<bool> BeforeChiR1;

	/* Parameter BeforeThetaR2 represents the value of the state
	* before Theta of round 2 (after Chi). */
	vector<bool> BeforeThetaR2;

	/* Vector AfterThetaR2 represents the value of the state
	* after Theta of round 2. */
	vector<bool> AfterThetaR2;

	/* Parameter BeforeChiR2 represents the value of the state
	* before Chi of round 2 (after Pi). */
	vector<bool> BeforeChiR2;

	/* Parameter BeforeThetaR3 represents the value of the state
	* before Theta of round 3 (after Chi). */
	vector<bool> BeforeThetaR3;

	/* Vector AfterThetaR3 represents the value of the state
	* after Theta of round 3. */
	vector<bool> AfterThetaR3;

	/* Parameter BeforeChiR3 represents the value of the state
	* before Chi of round 3 (after Pi). */
	vector<bool> BeforeChiR3;
};