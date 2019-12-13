#ifndef __HLBGrid_H__
#define __HLBGrid_H__

#include <vector>
#include "d:\\Program Files (x86)\\NAG\\FL26\\fldll26del\\c_headers\\nagmk26.h"
//#include "d:\\Program Files (x86)\\NAG\\FL25\\fldll254ml\\c_headers\\nagmk25.h"
//#include "d:\\Program Files (x86)\\NAG\\FL23\\flw3223dcl\\c_headers\\nagmk23.h"


namespace Grid
{
	int GetNumRows();
	int GetNumCols();

	//Getting and setting numbers of HLB
	double GetNumG(int irow, int jcol, int Health); //get total number of specified type 0 healthy and 1 infected
	double GetResNumG(int irow, int jcol, int Res); //get total number of resistence type 'res;
	double GetResSexNumG(int irow, int jcol, int Res, int Sex); //get total number of resistence type 'res' and Sex type 'Sex'
	double GetHealthNumG(int irow, int jcol,int Heal);  //get total number of resistence type 'heal'
	double GetSexNumG(int irow, int jcol, int Sex); //get total number of resistence type 'sex'
	void SetNumG(int irow, int jcol, int health, double NewNum); //set total number of bugs
	//int GetNumResTypes(); //Get number of resistant types. Currently set to 3 
	void SetTempBug(int irow, int jcol, int Health, double NewNum);
	double GetTempBug(int irow, int jcol, int Health);
	
	//Keeps a vector of larvae MAY NOT NEED
	void SetOWLav(int irow, int jcol, double totHPop, double totIPop); //get numbers of Bugs after they have grown and moved
	void GetOWLav(int irow, int jcol, int iyear, double& totHPop, double& totIPop); //get numbers of Bugs
	
	
	//Seting trees status
	void SetCrop(int irow, int jcol, int myStat, double per); // sets proportion of trees at a given status
	double GetCrop(int irow, int jcol, int myStat);
	int GetMaxNumStatus(); // retuns the number of status a tree can be in , i.e. health, infected, sympomatic, removed
	
	//These are related to number of cells in plantations might need revising for our needs
	int GetPlantaInfoSize(); //Returns the number of plantations in the landscape
	int GetNumCells(int PlantaCount); //returns number of grid cells in a plantation. Probably always 4
	int GetCellLocation(int PlantaCount, int cellCount); // thsi returns a number representing location icol*numRows+irow
	
	//Set up landscape comands
//	void SetLand(int newL, double propBelief, int& iflag); //if newL=1 then we have a set template Landscale  otherwise read in from saved file. PropBelief the proportion that believe in control
	void SetLand(int newL, double alpha_W, double beta_W, double alpha_B, double beta_B, int& iflag); //double alpha_W, double beta_W, double alpha_B, double beta_B are paraters for beta dristribution describing initial belief structures 

	
	
	void SetPlanta(int irow, int jcol, int Ctype);
	void SetBlock(int irow, int jcol, int Ctype);
	void SetCHMA(int irow, int jcol, int Ctype);

	//Get landscape info
	int GetPlanta(int irow, int jcol);
	int GetCHMA(int irow, int jcol);
	int GetBlock(int irow, int jcol);
	int GetBlockCHMA(int BlockID); //Returns Chima number for a given block. Might need revising
	int GetBlockPlanta(int BlockID); //returns the Planta that the Block is in

	//Apply control or not
	int GetControl(int irow, int jcol, int mnthago); //Get control 0 = no control 1 = control, mnthago how many months ago we are looking
	void SetControl(int irow, int jcol, int Control); //Set control 0 = no control 1 = control

	
	void ClearData();
	void ClearTempBugs();
	
	void SetIntVars(double jdist, double gridLen, double lambda);
	void GetIntVars(double& jdist, double& gridLen, double& lambda);
	//int GetNYKeep();
		
	//Comands to track invasion may not need
	void SetInvasion(int irow, int jcol, int yr);
	int GetInvasion(int irow, int jcol);
	
	//Neighbour info may not need all this
	void SetNeighbourlist(int myBlockNum, int NewBlockID); //creates a list of Block neighbours
	void GetNeighbourlist(int myBlockID, int myNeighbs[], int& NumNeigbs); //gets list of Block neighbours
	void DeleteNeighbourlist(int myBlockID);
	void SetConnectType(int myBlockID, int myCon);
	int GetConnectType(int myBlockID);

	void GetPlantaScore(int Plantcount, double& Score, double& BeliefInControl); //Gets plantation score or worry factor
	void SetPlantaScore(int Plantcount, double Score);  //Sets plantation score or worry factor
	void SetPlantaControlBelief(int Plantcount, double myBelief);  //Sets plantation score or worry factor
	void GetPlantaDisStatus(int Planta, int& mystat); //Stays whether sympotoms have gone over 5% in any cell in that plantation
	void SetPlantaDisStatus(int Planta, int mystat);
	void AddHealthLoss(int Plantcount, double incrInf); // cummulate the health loss across the plantation
	void ResetHealthLoss(int Plantacount);
	double GetHealthLoss(int Plantacount);


	void SetTotalPlantasFilled(int mycount);
	int GetTotalPlantasFilled();
	
}

namespace FlightP
{
	void ClearData();
	void Initialise();
	void GetBugFlight(std::vector<double>& MyDist);


}

double Integrate(int idist, int jdist, double lambda);


void GridIDToCoord(int Floc, int& irow, int& icol);
int CoordToGridID(int irow, int icol);

//void 
extern "C" void __stdcall fFunc2(const int& NDIM, const double x[], const int& NumFunc, double funcY[] );


#endif