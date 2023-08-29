#pragma once

#include "matrix_t.h"

#include <ilcplex/ilocplex.h>
#include <string> 
#include <iostream>
#include <ios>
#include <iomanip>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sstream>
#include <vector>
#include <algorithm>
#include <time.h>
#include <ctime>
#include <chrono>
#include <ctype.h>
#include <assert.h>
#include <limits>

class ProductionPlan
{
public:
	ProductionPlan(int, int, int, int, int, int); // Added integer for the Number of Facilities and Number of Stockpiles
	~ProductionPlan();

	void ReadData();
	void AssignBlockIDS();
	void AssignData(double, double, double, double, double, double, double*, double*, double*); //Updated according to accept the vector format of Processing Costs
	void CalculateBlockValues();
	void CreateModel();
	void AnalysisOfResults();

protected:

	MATRIX<int> mint;
	MATRIX<double> mdouble;
	MATRIX<char> mchar;

	double *bX, *bY, *bZ, *bTons, *bGrade, *blockValue, *wasteValue;
	int *blockID, *fromBlock, *toBlock;
	char *debugfilename;
	ofstream debugFile;

private:
	int nPeriods, nBlocks, nArcs, uIndicator, nFacilities, nStockpiles; // Added nFacilities and nStockpiles
	double metalPrice, refineCost, mining_cost, discount_rate, mining_cap_low, mining_cap_up;
	double  *processing_cap_low, *processing_cap_up, *processing_cost;
	double **processRecovery; // Added processRecovery to a two-dimensional array
	//Stockpile related member variables
	double **stockpile_recovery; // Added stockpileRecovery to a two-dimensional vector- pointer to pointer
	double *stockpile_capacity; // Added stockpileCapacity to a one-dimensional vector- pointer to pointer
};

