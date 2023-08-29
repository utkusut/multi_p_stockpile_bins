#include "ProductionPlan.h"

ProductionPlan::ProductionPlan(int periods, int blocks, int arcs, int units, int facilities, int stockpiles)
{
	system("copy 05_debugger.txt predebugger.txt");
	remove("debugger.txt");

	debugFile.open("05_debugger.txt", ios::app);
	debugFile.setf(ios::fixed, ios::floatfield);
	debugFile.precision(2);

	nPeriods = periods;
	nBlocks = blocks;
	nArcs = arcs;
	uIndicator = units;
	nFacilities = facilities;
	nStockpiles = stockpiles;

	metalPrice = refineCost = mining_cost = discount_rate = mining_cap_low = mining_cap_up = 0.0 ;
	processing_cost = nullptr;
    stockpile_capacity = nullptr;
    processing_cap_low = nullptr;
    processing_cap_up = nullptr;
    processRecovery = nullptr;
    stockpile_recovery = nullptr;

	mdouble.newMTX(bX, nBlocks);
	mdouble.newMTX(bY, nBlocks);
	mdouble.newMTX(bZ, nBlocks);
	mdouble.newMTX(bTons, nBlocks);
	mdouble.newMTX(bGrade, nBlocks);
	mdouble.newMTX(blockValue, nBlocks);
	mdouble.newMTX(wasteValue, nBlocks);
	mdouble.newMTX(processing_cost, nFacilities); // Initialize processing_cost with the correct dimensions
	mdouble.newMTX(processRecovery, nBlocks, nFacilities); // Initialize processRecovery with the correct dimensions
	mdouble.newMTX(stockpile_capacity, nStockpiles); // Initialize stockpile_capacity with the correct dimensions
	mdouble.newMTX(stockpile_recovery, nStockpiles, nFacilities); // Initialize stockpile_recovery with the correct dimensions
	mdouble.newMTX(processing_cap_low, nFacilities); // Initialize processing_cap_low with the correct dimensions
	mdouble.newMTX(processing_cap_up, nFacilities); // Initialize processing_cap_up with the correct dimensions
	

	mint.newMTX(blockID, nBlocks);
	mint.newMTX(fromBlock, nArcs);
	mint.newMTX(toBlock, nArcs);

}

ProductionPlan::~ProductionPlan(void)
{
	mdouble.delMTX(bX, nBlocks);
	mdouble.delMTX(bY, nBlocks);
	mdouble.delMTX(bZ, nBlocks);
	mdouble.delMTX(bTons, nBlocks);
	mdouble.delMTX(bGrade, nBlocks);
	mdouble.delMTX(blockValue, nBlocks);
	mdouble.delMTX(wasteValue, nBlocks);
	mdouble.delMTX(processing_cost, nFacilities); // Delete processingCosts with the correct dimensions
	mdouble.delMTX(processRecovery, nBlocks, nFacilities); // Delete processRecovery with the correct dimensions
	mdouble.delMTX(stockpile_capacity, nStockpiles); // Delete stockpile_capacity with the correct dimensions
	mdouble.delMTX(stockpile_recovery, nStockpiles, nFacilities); // Delete stockpile_recovery with the correct dimensions
	mdouble.delMTX(processing_cap_low, nFacilities); // Delete processing_cap_low with the correct dimensions
	mdouble.delMTX(processing_cap_up, nFacilities); // Delete processing_cap_up with the correct dimensions

	mint.delMTX(blockID, nBlocks);
	mint.delMTX(fromBlock, nArcs);
	mint.delMTX(toBlock, nArcs);
		debugFile.close();
}

void ProductionPlan::ReadData()
{
	debugFile << "Starting to read XYZ and block grades from 01_blockmodel.txt..."
		<< endl << endl << endl;

	std::cout << "Starting to read XYZ and block grades from 01_blockmodel.txt..."
		<< endl << endl << endl;

	ifstream myInput("01_blockmodel.txt", ios::in);
	ifstream myPrecedence("02_precedence.txt", ios::in);

	if (!myInput.is_open())
	{
		std::cerr << "Error: Unable to open file 01_blockmodel.txt for reading" << endl;
		return;
	}
		for (int i = 0; i < nBlocks; i++)
		{
			myInput >> bX[i] >> bY[i] >> bZ[i] >> bGrade[i] >> bTons[i];
			for (int p = 0; p < nFacilities; p++) // Read the processing recovery for each block and each processing facility
            {
                myInput >> processRecovery[i][p]; // where i is the block and p is the corresponding Processing Facility
            }
			//std::cout << (i + 1) << endl;
		}
}

	myInput.close();

	if (!myPrecedence.is_open())
	{
		std::cerr << "Error: Unable to open 02_precedence.txt for reading." << endl;
        return;
    }
		for (int i = 0; i < nArcs; i++)
		{
			myPrecedence >> fromBlock[i] >> toBlock[i];
			//std::cout << (i + 1) << endl;
		}
		myPrecedence.close();

		debugFile << "Finished reading XYZ and block grades from 01_blockmodel.txt."
		<< endl << endl << endl;

		std::cout << "Finished reading XYZ and block grades from 01_blockmodel.txt."
		<< endl << endl << endl;
}

void ProductionPlan::AssignBlockIDS()
{
	debugFile << "Assign IDs ..."
		<< endl << endl << endl;

	std::cout << "Assign IDs ....."
		<< endl << endl << endl;

	int k = 0;

	for (int j = 0; j < nBlocks; j++)
	{
		blockID[k] = k + 1;
		k++;
	}

	ofstream myOutput("03_InitialData.txt", ios::out);
	myOutput.setf(ios::fixed, ios::floatfield);

	    for (int i = 0; i < nBlocks; i++)
    {
        myOutput << setiosflags(ios::left) << setprecision(2)
            << setw(10) << blockID[i]
            << setw(10) << bX[i]
            << setw(10) << bY[i]
            << setw(10) << bZ[i]
            << setw(10) << bTons[i]
            << setw(10) << bGrade[i];
        
        for (int p = 0; p < nFacilities; p++)
        {
            myOutput << setw(10) << processRecovery[i][p];
        }

        myOutput << endl;
    }

    for (int i = 0; i < nArcs; i++)
    {
        myOutput << setiosflags(ios::left)
            << setw(10) << fromBlock[i]
            << setw(10) << toBlock[i] << endl;
    }
	myOutput.close(); 
	debugFile << "Assign IDs ...."
		<< endl << endl << endl;

	std::cout << "Assign IDs ....."
		<< endl << endl << endl;
}

void ProductionPlan::AssignData(double mP, double rC, double mC, double dR, double mLC, double mUC, double *pLC, double *pUC, double* pC)
{
	metalPrice = mP;
	refineCost = rC;
	mining_cost = mC;
	discount_rate = dR;
	mining_cap_low = mLC;
	mining_cap_up = mUC;


	//assign processing costs
	for (int p = 0; p < nFacilities; p++)
		{
		processing_cost[p] = pC[p];
		processing_cap_low[p]= pLC[p];
		processing_cap_up[p] = pUC[p];
		}
}

void ProductionPlan::CalculateBlockValues() //This code calculates the maximum block value across all processing facilities and stores the highest value for each block, taking into account the best processing facility for maximizing the block value.
{
	for (int i = 0; i < nBlocks; i++)
	{
		double maxValue = -std::numeric_limits<double>::max(); //This is the value of the block if it is treated as ore
		double wasteValue = (-1)*mining_cost*bTons[i]; //This is the value of the block if it is treated as waste

		for (int p = 0; p < nFacilities; p++) // Add the facility loop here
		{
			double value = 0.0;//This is the value of the block if it is treated as ore and processed at facility p
			if (uIndicator == 1)
			{
				value = ((metalPrice - refineCost) * (bGrade[i] / 100) * processRecovery[i][p] * bTons[i]) - (mining_cost * bTons[i]) - (processing_cost[p] * bTons[i]);
			}
		
			else if (uIndicator == 2)
			{
				value = ((metalPrice - refineCost) * bGrade[i] * processRecovery[i][p] * bTons[i]) - (mining_cost * bTons[i]) - (processing_cost[p] * bTons[i]);
			}
			if (value > maxValue)
			{
				maxValue = value;
			}
		}
			if (maxValue <= 0.0)
			{
				blockValue[i] = wasteValue;
			}
			else
			{
				blockValue[i] = maxValue;
			}	
	
		
		debugFile << setiosflags(ios::left) << setw(10) << (i + 1) << blockValue[i] << endl;
	}

	debugFile << endl;
}

void ProductionPlan::CreateModel()
{
	IloEnv env;
	IloModel mod(env);
	IloTimer solutiontime(env);

	typedef IloArray<IloArray<IloNumVarArray> > Array3D;

	Array3D x(env, nBlocks); // Decision variable for the production of each block at each facility in each period

	for(int i = 0; i < nBlocks; i++) 
	{
		x[i] = IloArray<IloNumVarArray>(env, nFacilities); 

		for (int p = 0; p < nFacilities; p++)
		{

		x[i][p] = IloNumVarArray(env, nPeriods, 0, 1, ILOINT); 

		}
	}
	
	for (int i = 0; i < nBlocks; i++)
	{
		for (int p = 0; p < nFacilities; p++)
		{
			for (int t = 0; t < nPeriods; t++) 

			debugFile << i << "  " << p << " " << t << "  " << x[i][p][t] << endl; 
		}
	}
	
	typedef IloArray<IloNumVarArray> Array2D;
	
	Array2D W(env, nBlocks); // Decision variable for treating each block as waste in each period

		for(int i = 0; i < nBlocks; i++) 
		{
			W[i] = IloNumVarArray(env, nPeriods, 0, 1, ILOINT); 
		}


	debugFile << endl << endl;

IloExpr objExpression(env); //Objective Function Expression

	int n = 0; //Loop integer, Loops until the objective completed.

	// Objective Function - Discounted value of the production plan (maximize) 
	for (int i = 0; i < nBlocks; i++) {
		for (int p = 0; p < nFacilities; p++) {
			for (int t = 0; t < nPeriods; t++) {
				n++;
				objExpression += (blockValue[i] * x[i][p][t]) / (pow((1 + (discount_rate / 100)), (t + 1)));
				debugFile << n << "  " << i << "  " << p << "  " << t << "  " << x[i][p][t] << endl;
			}
		}		
	}

	for (int i = 0; i < nBlocks; i++) {
		for (int t = 0; t < nPeriods; t++) {
			n++;
			objExpression += (wasteValue * W[i][t]) / (pow((1 + (discount_rate / 100)), (t + 1)));
			debugFile << n << "  " << i << "  WASTE  " << t << "  " << W[i][t] << endl;
		}
	}

		debugFile << endl << endl;

IloObjective objFunction(env, objExpression, IloObjective::Maximize);

mod.add(objFunction);

	
		for (int i = 0; i < nBlocks; i++) 
	{
		IloExpr reserveExpression(env);

		for (int t = 0; t < nPeriods; t++) 
		{
			for (int p = 0; p < nFacilities; p++) 
			{
				reserveExpression += x[i][p][t];
			}
			reserveExpression += W[i][t];
		}

		mod.add(reserveExpression <= 1); // Reserve Constraint - Each Block can be mined at once!
	}

	 //Processing Facility Selection Constraint - Each block can assign only one processing facility
	for (int i = 0; i < nBlocks; i++)
	{
		for (int t = 0; t < nPeriods; t++)
		{
			IloExpr facilitySelection(env);

			for (int p = 0; p < nFacilities; p++)
			{
				facilitySelection += x[i][p][t];
			}

			mod.add(facilitySelection = 1);
		}
	}
	

	for (int t = 0; t < nPeriods; t++)
{
    IloExpr miningExpression(env);

    for (int i = 0; i < nBlocks; i++)
    {
        for (int p = 0; p < nFacilities; p++) 
        {
            miningExpression += (bTons[i] * x[i][p][t]); 
        }
        miningExpression += (bTons[i] * W[i][t]); // Adding the tons treated as waste
    }

    mod.add(miningExpression >= mining_cap_low); // Lower Mining Capacity Constraint
    mod.add(miningExpression <= mining_cap_up); // Upper Mining Capacity Constraint
}

	
	for (int t = 0; t < nPeriods; t++)
	{
		for (int p = 0; p < nFacilities; p++)
		{
			IloExpr processExpression(env);

        for (int i = 0; i < nBlocks; i++)
        {
            if (blockValue[i] > 0)
            {
                processExpression += (bTons[i] * x[i][p][t]);
            }
        }

        mod.add(processExpression >= processing_cap_low[p]);
		mod.add(processExpression <= processing_cap_up[p]);
		}
	}
	   
		int fBlock = 0, tBlock = 0;

		for (int j = 0; j < nArcs; j++)
		{
			fBlock = fromBlock[j]; // This is the preceding block b'
			tBlock = toBlock[j];   // This is the current block b

			for (int t = 0; t < nPeriods; t++)
			{
				for (int p = 0; p < nFacilities; p++) // Loop over facilities
				{
					IloExpr sumExpression(env);

					// Sum over all preceding blocks and all facilities for the current block
					for (int tPrime = 0; tPrime <= t; tPrime++)
					{
						sumExpression += x[fBlock - 1][p][tPrime];
						sumExpression += W[fBlock - 1][tPrime];
					}

					// Add the constraint for the block going to processing facility
					mod.add(n * x[tBlock - 1][p][t] - sumExpression <= 0);

					// Add the constraint for the block going to waste
					mod.add(n * W[tBlock - 1][t] - sumExpression <= 0);
				}
			}
		}


		debugFile << endl;
	}

	

	IloCplex cplex(env);
	cplex.extract(mod);
	cplex.exportModel("04_Formulation.lp");
	
	cplex.setParam(IloCplex::ClockType, 1);
	cplex.setParam(IloCplex::EpGap, 0.01);
	cplex.setParam(IloCplex::TiLim, 2592000);

	auto start = std::chrono::system_clock::now();
	solutiontime.start();
	cplex.solve();
	solutiontime.stop();
	auto end = std::chrono::system_clock::now();

	std::chrono::duration<double> elapsed_seconds = end - start;

	ofstream myOutput("05_OptimalPP.txt", ios::out);
	myOutput.setf(ios::fixed, ios::floatfield);

	if (myOutput.is_open())
	{
		IloNum objValue = cplex.getObjValue();

		myOutput << "Discounted value of the production plan = " << setprecision(0) << objValue << endl;
		myOutput << "Solution time = " << setprecision(5) << solutiontime.getTime() << " seconds" << endl << endl;
		myOutput << "Solution time = " << setprecision(5) << elapsed_seconds.count() << " seconds" << endl << endl;

		for (int t = 0; t < nPeriods; t++)
		{
			for (int i = 0; i < nBlocks; i++)
			{
				for (int p = 0; p < nFacilities; p++) // Added this loop for processing facilities

					if (cplex.getValue(x[i][p][t]) > 0)
					{
						myOutput << setiosflags(ios::left) << setprecision(0)
							<< setw(10) << bX[i]
							<< setw(10) << bY[i]
							<< setw(10) << bZ[i]
							<< setw(10) << bTons[i]
							<< setprecision(5)
							<< setw(15) << bGrade[i]
							<< setw(10) << setprecision(0) << (t+1)
							<< setw(10) << p + 1 // Added facility number (p + 1) to the output
							<< setw(10) << cplex.getValue(x[i][p][t]) << endl;
				}
			}
		}

		ofstream myOutput1("06_PitProcess.txt", ios::out); // Print out the blocks that are in the UPL and goes processing
		myOutput1.setf(ios::fixed, ios::floatfield);

		for (int i = 0; i < nBlocks; i++)
		{
			for (int p = 0; p < nFacilities; p++) // Added this loop for processing facilities
			{
				for (int t = 0; t < nPeriods; t++)
				{
					if (cplex.getValue(x[i][p][t]) > 0)
					{
						myOutput1 << setiosflags(ios::left) << setprecision(0)
							<< setw(10) << blockID[i]
							<< setw(10) << (p)
							<< setw(10) << (t) << endl;
					}
				}
			}			
		}
	}
}

void ProductionPlan::AnalysisOfResults()
{
	ofstream resultsFile("07_results.txt", ios::out);
	resultsFile.setf(ios::fixed, ios::floatfield);

	resultsFile << setiosflags(ios::left) << setprecision(0);

	ofstream modelFile("08_blockmodelSchedule.txt", ios::out);
	modelFile.setf(ios::fixed, ios::floatfield);

	ofstream blockToFacilityFile("09_blockToFacility.txt", ios::out);
	blockToFacilityFile.setf(ios::fixed, ios::floatfield);

	ofstream processedMaterialFile("10_processedMaterial.txt", ios::out);
	processedMaterialFile.setf(ios::fixed, ios::floatfield);

	//Pit to Process
	ifstream pToP("06_PitProcess.txt", ios::in);

	int *mbID, *mTime, *mFacility;

	mint.newMTX(mbID, (nBlocks * nPeriods));
	mint.newMTX(mFacility, (nBlocks * nFacilities));
	mint.newMTX(mTime, (nBlocks * nPeriods));
	int i = 0, nPitProcess = 0;
	while (pToP >> mbID[i] >> mFacility[i] >> mTime[i])
	{
		i++;
	}

	nPitProcess = i;

	for (int b = 0; b < nBlocks; b++)
	{
		for (int i = 0; i < nPitProcess; i++)
		{
			if (blockID[b] == mbID[i])
			{
				modelFile << setiosflags(ios::left) << setprecision(0)
					<< setw(10) << bX[b]
					<< setw(10) << bY[b]
					<< setw(10) << bZ[b]
					<< setw(10) << bTons[b]
					<< setprecision(5)
					<< setw(15) << bGrade[b]
					<< setw(15) << (mFacility[i] + 1) << endl
					<< setw(15) << (mTime[i] + 1) << endl;

				goto out_block;
			}
		}
	
		modelFile	<< setiosflags(ios::left) << setprecision(0)
					<< setw(10) << bX[b]
					<< setw(10) << bY[b]
					<< setw(10) << bZ[b]
					<< setw(10) << bTons[b]
					<< setprecision(5)
					<< setw(15) << bGrade[b]
					<< setw(15) << (nPeriods + 100) << endl;

		out_block: debugFile << endl;
	}

	double *qtWaste, *qtOre, *qtMetal, *cutG, *avgG, *cashF, *netPV, *totalProcessed;

	mdouble.newMTX(qtWaste, nPeriods);
	mdouble.newMTX(qtOre, nPeriods);
	mdouble.newMTX(qtMetal, nPeriods);
	mdouble.newMTX(cutG, nPeriods);
	mdouble.newMTX(avgG, nPeriods);
	mdouble.newMTX(cashF, nPeriods);
	mdouble.newMTX(netPV, nPeriods);
	mdouble.newMTX(totalProcessed, nFacilities);

	for (int p = 0; p < nFacilities; p++)
	{
		totalProcessed[p] = 0.0;
	}

	for (int i = 0; i < nPitProcess; i++) // loop over pit process entries
	{
		int b = mbID[i];
		int p = mFacility[i];
		int t = mTime[i];

		totalProcessed[p] += bTons[b];

		blockToFacilityFile << setiosflags(ios::left) << setprecision(0)
			<< setw(10) << "Block " << blockID[b] << "is sent to Processing Facility "
			<< setw(10) << (p + 1)<< "in period "
			<< setw(10) << (t + 1) << "." << endl;
	}

	for (int p = 0; p < nFacilities; p++)
	{
		processedMaterialFile << setiosflags(ios::left) << setprecision(0)
			<< setw(10) << "Processing Facility " << (p + 1) << " processes "
			<< setw(10) << totalProcessed[p] << " tons of material" << endl;
	}

	for (int t = 0; t < nPeriods; t++)
	{
			qtWaste[t] = 0.0;
			qtOre[t] = 0.0;
			qtMetal[t] = 0.0;
			cutG[t] = 1000000000.0;
			avgG[t] = 0.0;	
			netPV[t] = 0.0;
	}

	for (int i=0: i<nPitProcess; i++) // Instead of the below code looping over all blocks,periods,facilities; loop over pit process entries
	{
		double sumG = 0.0; 
		int b = mbID[i];
		int p = mFacility[i];
		int t = mTime[i];

		if (blockValue[b] <= 0)
		{
			qtWaste[t] += bTons[b];
		}
		else
		{
			qtOre[t] += bTons[b];

			for (int p = 0; p < nFacilities; p++)
			{
				qtMetal[t] += (processRecovery[b][p] * bTons[b]);
			}
			sumG += (bGrade[b] * bTons[b]);

			if (bGrade[b] < cutG[t])
			{
				cutG[t] = bGrade[b];
			}
		}
	}
	/*for (int t = 0; t < nPeriods; t++)
	{
		double sumG = 0.0;

		for (int b = 0; b < nBlocks; b++)
		{
			for (int i = 0; i < nPitProcess; i++)
			{
				if (blockID[b] == mbID[i] && t == mTime[i])
				{
					if (blockValue[b] <= 0)
					{
						qtWaste[t] += bTons[b];
					}
					else
					{	
						qtOre[t] += bTons[b];

						for (int p = 0; p < nFacilities; p++)
						{
							qtMetal[t] += (processRecovery[b][p] * bTons[b]);
						}
						sumG += (bGrade[b] * bTons[b]);

						if (bGrade[b] < cutG[t])
						{
							cutG[t] = bGrade[b];
						}
					}
				}
			}
		}*/
	for (int t = 0; t < nPeriods; t++)
	{	
		avgG[t] = (sumG / qtOre[t]);

		if (uIndicator == 1)
		{
			qtMetal[t] = (qtMetal[t] * (avgG[t] / 100));
		}

		if (uIndicator == 2)
		{
			qtMetal[t] = (qtMetal[t] * avgG[t]);
		}

		double totalProcessCost = 0.0;
		for (int i = 0; i < nPitProcess; i++) // loop over pit process entries
		{
			int b = mbID[i];
			int p = mFacility[i];

			if(t == mTime[i])
			{
				totalProcessCost += (processing_cost[b][p] * bTons[b]);
			}
		}

		cashF[t] = ((metalPrice - refineCost) * qtMetal[t]) - mininig_cost * (qtWaste[t] + qtOre[t]) - totalProcessCost;
	}

	for (int t = 0; t < nPeriods; t++)
	{
		for (int t1 = t; t1 < nPeriods; t1++)
		{
			if (t == 0)
			{
				netPV[t] += (cashF[t1] / (pow((1 + (discount_rate / 100)), (t1 + 1))));
			}
			else
			{
				netPV[t] += (cashF[t1] / (pow((1 + (discount_rate / 100)), (t1 - t + 1))));
			}			
		}
	}

	resultsFile << endl;

	resultsFile << setw(10) << "Year" 
			    << setw(20) << "Quantity of Waste" 
				<< setw(20) << "Cut-off Grade"
				<< setw(20) << "Average Grade"
				<< setw(20) << "Quantity of Ore" 
				<< setw(20) << "Quantity of Metal" 
				<< setw(20) << "Cash Flow" 
				<< setw(20) << "NPV" << endl;

	resultsFile << setw(10) << ""
				<< setw(20) << "(tonnes)"
				<< setw(20) << "(%)"
				<< setw(20) << "(%)"
				<< setw(20) << "(tonnes)" 
				<< setw(20) << "(tonnes)" 
				<< setw(20) << "($)" 
				<< setw(20) << "($)" << endl << endl;

	for (int t = 0; t < nPeriods; t++)
	{
		resultsFile << setw(10) << (t+1)
					<< setw(20) << qtWaste[t]
					<< setw(20) << setprecision(5) << cutG[t]
					<< setw(20) << setprecision(5) << avgG[t]
					<< setw(20) << setprecision(0) << qtOre[t] 
					<< setw(20) << setprecision(0) << qtMetal[t] 
					<< setw(20) << setprecision(0) << cashF[t] 
					<< setw(20) << setprecision(0) << netPV[t] << endl;
	}

	mint.delMTX(mbID, (nBlocks * nPeriods));
	mint.delMTX(mTime, (nBlocks * nPeriods));
	mint.delMTX(mFacility, (nBlocks * nPeriods));

	mdouble.delMTX(qtWaste, nPeriods);
	mdouble.delMTX(qtOre, nPeriods);
	mdouble.delMTX(qtMetal, nPeriods);
	mdouble.delMTX(cutG, nPeriods);
	mdouble.delMTX(avgG, nPeriods);
	mdouble.delMTX(cashF, nPeriods);
	mdouble.delMTX(totalProcessed, nFacilities);
	
}
