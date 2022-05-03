#include <iostream>
#include <stdio.h>
#include <fstream>
#include <algorithm> 
#include <chrono> 
#include <stdlib.h>
#include <malloc.h>
#include "Subterranean.h"
#include "TrailSeedIterator.h"

using namespace std;
/* To measure the compiling time. */
using namespace std::chrono;
int main(int argc, char* argv[]) {
	auto start = high_resolution_clock::now();
	/* A machine readable file of required 3-round trail cores. */
	ofstream fout3Machine;
	/* A file of required 3-round trail cores. */
	ofstream fout3;
	/* A file of the histogram of required 3-round trail cores. */
	ofstream foutHisto3;
	/* A file of all 4-round trail cores up to the given certain weight. */
	ofstream fout4;
	/* A file of the histogram of all 2-round trail cores up to the given certain weight. */
	ofstream foutHisto2;
	/* A file of the histogram of all 2-round trail cores up to the given certain weight in details. */
	ofstream foutHisto2Det;
	fout3Machine.open("3-r.csv");
	fout3.open("3-rTrailCores.txt");
	fout4.open("4-rTrailCores.txt");
	foutHisto3.open("Histo3-r.txt");
	foutHisto2.open("Histo2-r.txt");
	foutHisto2Det.open("Histo2-r_det.txt");


	/* Maximum weight of the 4-round trial cores. */
	// It will take almost 2 weeks to scan the space up to 57.
	int maxWeightFourRound = 57;
	/* Maximum weight of the 2-round trial cores. */
	int maxWeightTwoRound = maxWeightFourRound / 2;
	/* Array to store the histogram of all 2-round trail cores up to the certain weight (maxWeightTwoRound). */
	int his2Rout[100] = {};
	/* Array to store the histogram of some 3-round trail cores. */
	int his3Rout[100] = {};
	/* Array to store the histogram of all 4-round trail cores up to the certain weight (maxWeightFourRound). */
	int his4Rout[100] = {};

	/* We first backward generate all 2-round trail cores and extend them in the backward direction,
	* and then we generate all 2-round trail cores and extend them in the forward direction. */
	for (int aForwardBackward = 0; aForwardBackward <= 1; aForwardBackward++) { 
		/** TrailSeedIterator constructor.
		* The first parameter represents the direction of extension,
		* the second parameter represents the maximum weight of 2-round trail cores, 
		* the third parameter represents the maximum weight of 4-round trail cores, 
		* the remaining parameters represent the arrays for storing the histograms. */
		TrailSeedIterator T(aForwardBackward, maxWeightTwoRound, maxWeightFourRound, his2Rout, his3Rout, his4Rout);
		/* Tree search.*/
		while (T.next()) {
			/* Method to generate required 3-round trail cores.*/
			T.Three_RoundTrailCore();
		}
		/* Method to show required 3-round trail cores in the output.*/
		T.showw3(fout3);
		/* Method to generate all 4-round trail cores up to maxWeightFourRound.*/
		T.Four_RoundTrailCore();
		/* Method to show required 4-round trail cores in the output.*/
		T.showw4(fout4);

		/* Copy the histograms of backward direction
		* in order to add the histograms of the forward direction to it.*/
		for (int i = 0; i < 100; i++) {
			his2Rout[i] = T.his2R[i];
			his3Rout[i] = T.his3R[i];
			his4Rout[i] = T.his4R[i];
		}
		/* Generate histograms if both forward and backward extension is done.*/
		if (aForwardBackward == 1) {
			T.hi(foutHisto3);
			T.hi2(foutHisto2);
			T.hi22(foutHisto2Det);
		}
		/* Generate a machine readable file of required 3-round trail cores.*/
		T.machinReadable(fout3Machine);
	}
	/* Computing the compiling time.*/
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<hours>(stop - start);
	fout4 << "\n\nExecution time: " << duration.count() << " hours" << endl;

	fout3.close();
	fout4.close();
	foutHisto3.close();
	foutHisto2.close();
	foutHisto2Det.close();
	fout3Machine.close();
	return 0;
}