#include "HLBGrid.h"
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <direct.h>




int const NumCropTypes=4;
int const MaxNumNeighbours=50;
int const MaxNumPCells=12169; //maximum number cells per plantation

typedef double ResHealthSex[2];// index: 0=healthy, 1=infected
								

typedef double ThreeDoubles[3];

class PInfo
{
public:
	PInfo()
	{
	//	BlockNum=-1;
		DisStatus=0; //We assume the woory level is zero on initialisation
		Score=0;
		BeliefInControl = 0; 
		MonthlyLossofHealth = 0;
	//	numCells=0;
		//for (int icount=0; icount<MaxNumPCells; icount++)
		//{
		Cells.clear();
		//}

		for (int icount=0; icount<MaxNumNeighbours; icount++)
		{
			NeighbourBlockID[icount]=-1; // list of neighbours
		}
		NumNeighbours=0; //number of neighbours
		ConnectionType=-1;
		UniquePlantaID=-1;
		

	}
	PInfo(int myCHMA, int myPlantaID,  int mySym, double myscore, double mybelief, std::vector<int> mycells)
	{
		FCHMA=myCHMA;
	//	BlockNum=-1; //not used
		UniquePlantaID=myPlantaID;
		DisStatus=mySym; //0 if < 5% sympotmatic and 1 if more than that in any one cell
		Score=myscore;
		BeliefInControl = mybelief;
		MonthlyLossofHealth = 0;
		//numCells=myNumcells;
		Cells.clear();
		for (int icount=0; icount<mycells.size(); icount++)
		{
			Cells.push_back(mycells[icount]);
		}
		for (int icount=0; icount<MaxNumNeighbours; icount++)
		{
			NeighbourBlockID[icount]=-1; // list of neighbours
		}
		NumNeighbours=0; //number of neighbours
		ConnectionType=-1;
		
	}
	~PInfo()
	{
		NumNeighbours=0; //number of neighbours
		Cells.clear();
	}
	void SetBlockScores(double myScore) // this is planta not block 
	{
		Score=myScore;
	}
	void SetPBelief(double mybelief)
	{
		BeliefInControl = mybelief;
	}
	void SetNeighbour(int NID)
	{
		if (NumNeighbours<MaxNumNeighbours) 
		{
			bool alreadyCounted=false;
			for (int icount=0; icount<NumNeighbours; icount++)
			{
				if (NeighbourBlockID[icount]==NID)
						alreadyCounted=true;
			}
			if (!alreadyCounted)
			{
				NeighbourBlockID[NumNeighbours]=NID;
				NumNeighbours++;
			}
		}
		else
		{
			throw std::logic_error("too many neighbours");
		}
	}
	void SetConnect(int myCon)//This tells us how the Block is connected to the information highway. -1 = no connection, 0=neighbour; 1=Planta; 2=CHMA
	{
		ConnectionType=myCon;
	}
	
	int GetnumCells(){return Cells.size();}
	int GetCellID(int CellCount){return Cells[CellCount];}
	int GetnumNeighbours(){return NumNeighbours;}
	int GetNeighbourID(int ncount){return NeighbourBlockID[ncount];}
	double GetScore(){return Score;}
	double GetBelief(){ return BeliefInControl; }
	int GetFmCHMA(){return FCHMA;}
//	int GetFmPlanta(){return UniquePlantaID;}
	void GetNeighbour(int NID, int myNeighbs[], int& NumNeibs)
	{
		for (int icount=0; icount<NumNeighbours; icount++)
		{
			myNeighbs[icount]=NeighbourBlockID[icount];
		}

		NumNeibs=NumNeighbours;
		
	}
	int GetPlantNum(){return UniquePlantaID;}
	void DeleteNeigh(){NumNeighbours=0;}
	int GetConnect()//This tells us how the Block is connected to the information highway. -1 = no connection, 0=neighbour; 1=Planta; 2=CHMA
	{
		return ConnectionType;
	}
	int GetDisStatus(){ return DisStatus;}
	void SetDisStatus(int myStat){DisStatus=myStat;}
	void AddHLoss(double exloss){ MonthlyLossofHealth = MonthlyLossofHealth + exloss; }
	void ResetHealthLoss(){ MonthlyLossofHealth = 0; }
	double GetHLoss(){ return MonthlyLossofHealth; }
	

private:
	//	int BlockNum;
		int DisStatus; //If symptoms above 5% somewhere on the plantation then worry sets in
		double Score; //this is the 'worry' 0 is not worried about infection and 1 is very worried
		double BeliefInControl; //If a farmer firmly believes in the control strategy then we score him 1 is not zero. 
		std::vector<int> Cells; // field locations field=icols*maxnumrows+irows
		//int numCells; // number of cells in a Block
		int FCHMA;
		int NeighbourBlockID[MaxNumNeighbours]; // list of neighbours
		int NumNeighbours; //number of neighbours
		int ConnectionType; //This tells us how the Block is connected to the information highway. -1 = no connection, 0=neighbour; 1=Planta; 2=CHMA
        int UniquePlantaID; //the Planta the Block is in
		double MonthlyLossofHealth;
		
};

class TwoVars
{
public:
  TwoVars(){Xs[0]=0; Xs[1]=0;};
  TwoVars(double val1, double val2){Xs[0]=val1; Xs[1]=val2;};
  ~TwoVars(){};
  TwoVars(const TwoVars& ATheX){Xs[0]=ATheX.Xs[0]; Xs[1]=ATheX.Xs[1];};
  TwoVars& operator=(const  TwoVars& ATheX){Xs[0]=ATheX.Xs[0]; Xs[1]=ATheX.Xs[1]; return *this;};
  bool  operator<(const  TwoVars& ATheX)
  {
    bool A=false; 
    if(Xs[0]<ATheX.Xs[0]) A=true;
    return A;
  };
  bool  operator>(const  TwoVars& ATheX){bool A=false; if(Xs[0]>ATheX.Xs[0]) A=true; return A;};
  double Xs[2];
};


class Bugs
{
public:
	
	Bugs()
	{
			BugT[0]=0; //Healthy
			BugT[1]=0; //Sick
		
	}

	~Bugs()
	{
		TotalOWLavI.clear(); 
		TotalOWLavH.clear(); 
		BugT[0]=0; //Healthy
		BugT[1]=0; //Sick

	}
	double GetNum(int Health){return BugT[Health];} //get total number of specified type
	
	void SetNum(int Health, double NewNum){BugT[Health]=NewNum;} 
	
	void SetOWLav(double totHPop, double totIPop){TotalOWLavH.push_back(totHPop); TotalOWLavI.push_back(totIPop);}
	
	void GetOWLav(int iyear, double& totHPop, double& totIPop){totHPop=TotalOWLavH[iyear]; totIPop=TotalOWLavI[iyear];}
	
	void ClearTempBug()
	{
		TempBug[0]=0; //healthy
		TempBug[1]=0; //sick
		
	}
	void Clear()
	{
		TotalOWLavI.clear(); 
		TotalOWLavH.clear(); 
		BugT[0]=0; //Healthy
		BugT[1]=0; //Sick

	}
	void SetTempBug( int Health, double NewNum)
	{
		TempBug[Health]=NewNum;
	}
	double GetTempBug(int Health)
	{
		return TempBug[Health];
	}


private:
	ResHealthSex BugT; // number of each type of adult Bug (average per plant/unit area)
	std::vector<double> TotalOWLavH;
	std::vector<double> TotalOWLavI;
	ResHealthSex TempBug; //a structure to keep track of travelling Bugs
	
};


namespace Grid
{
	int const HistLen(12);
	int const Nrows = 350; //<-CHMA IR //290; //CHMA 15 // 200; // dummy //265;  //CHMA 54 //
	int const Ncols = 400; //<-CHMA IR //880; //CHMA 15 // 200; // dummy //220; //CHMA 54 //
	int const NStatus=5; //Status of trees. This is either Healthy=0, Infected but not infectious = 1, infected infectctious no symtoms=2, Symptoms=3, Removed (which includes the area never planted) =4 
	int const MaxNumYears(200);
	Bugs GridBugs[Nrows][Ncols];
	double Trees[Nrows][Ncols][NStatus]; //Holds proportion of trees in each CHMA listed under status
	int G_Control[Nrows][Ncols][HistLen]; //We now hold control for 12 months. the zero position is control for last month the 1 position for the month before
	//Note that  areas that cannot be cropped are 4 with fix=1

	double IntVars[3]; //keeps parameters for dispersal integration can be made bigger
	//each cell has a CHMA, Planta and Block ID. The combination is unique.
	int CHMAID[Nrows][Ncols];
	int PlantaID[Nrows][Ncols];
	int BlockID[Nrows][Ncols];
	int invasionYr[Nrows][Ncols]; // this is to store the year that the invasion of the Bug becomes > a fixed number
	std::vector<PInfo> PlantationInfo; //This stores the info about each Block. Uniquie Block ID is offset by one from indexing, i.e. Block 1 info is held in index 0
	int TotalPlantasFilled;

	int GetTotalPlantasFilled(){return TotalPlantasFilled;}

	int GetNumRows(){return Nrows;}
	
	int GetNumCols(){return Ncols;}

	//int GetNumResTypes(){return NumResTypes;}
	
    double GetNumG(int irow, int jcol, int Health)
	{
		if ((irow<GetNumRows())&&(jcol<GetNumCols())&&(irow>-1)&&(jcol>-1))
		{
			return GridBugs[irow][jcol].GetNum(Health);
		}
		else
			throw std::logic_error("Error in GetNumG");
	
	} //get total number of specified type
	
	int GetPlantaInfoSize(){return PlantationInfo.size();}

	int GetNumCells(int BlockID){return PlantationInfo[BlockID].GetnumCells();}

	int GetCellLocation(int PlantaID, int cellCount){ return PlantationInfo[PlantaID].GetCellID(cellCount);}

	void GetPlantaScore(int PlantaCount, double& myScore, double& myBelief) //Gets plantation score or worry factor
	{
		myScore=PlantationInfo[PlantaCount].GetScore();
		myBelief = PlantationInfo[PlantaCount].GetBelief();
	}

	void GetPlantaDisStatus(int PlantaCount, int& myStat)
	{
		if (PlantationInfo[PlantaCount].GetPlantNum() == PlantaCount + 1)
			myStat = PlantationInfo[PlantaCount].GetDisStatus();
		else
			throw std::logic_error("PlantaCount and index don't match GPDS"); //the function should not exit here unless there is no Maize
	}

	void SetNeighbourlist(int myBlockID, int NewBlockID) //creates a list of Block neighbours
	{
		if ((NewBlockID!=0)&&(NewBlockID!=myBlockID))
			PlantationInfo[myBlockID].SetNeighbour(NewBlockID);
	}

	void GetNeighbourlist(int myBlockID, int myNeighbs[], int& NumNeibs) //gets list of Block neighbours
	{
			PlantationInfo[myBlockID].GetNeighbour(myBlockID, myNeighbs, NumNeibs);
			
	}

		

	void DeleteNeighbourlist(int myBlockID)
	{
		PlantationInfo[myBlockID].DeleteNeigh();
	}


	void SetTotalPlantasFilled(int myCounts){ TotalPlantasFilled=myCounts;}

	void SetNumG(int irow, int jcol, int Health, double NewNum)
	{ 
		if ((irow<GetNumRows())&&(jcol<GetNumCols()))
		{
			GridBugs[irow][jcol].SetNum(Health, NewNum);
		}
		else
			throw std::logic_error("Error initialising in grid cells that don't exist");

	} 

	void SetConnectType(int myBlockID, int myCon)
	{
		PlantationInfo[myBlockID].SetConnect(myCon);
	}

	int GetConnectType(int myBlockID)
	{
		return PlantationInfo[myBlockID].GetConnect();
	}

	void SetOWLav(int irow, int jcol, double totHPop, double totIPop){ GridBugs[irow][jcol].SetOWLav(totHPop, totIPop);}
	
	
	

	void SetPlanta(int irow, int jcol, int Ctype)
	{
		if ((irow<GetNumRows())&&(jcol<GetNumCols())&&(irow>-1)&&(jcol>-1))
		{
			PlantaID[irow][jcol]=Ctype; 
		}
		else
			throw std::logic_error("Error in SetPlanta");
			//int junk=1;

	}


	


	void SetBlock(int irow, int jcol, int Ftype)
	{
		if ((irow<GetNumRows())&&(jcol<GetNumCols())&&(irow>-1)&&(jcol>-1))
		{
			BlockID[irow][jcol]=Ftype; 
		}
		else
			throw std::logic_error("Error in SetBlock");
			//int junk=1;

	}


	void SetCHMA(int irow, int jcol, int Stype)
	{
		if ((irow<GetNumRows())&&(jcol<GetNumCols())&&(irow>-1)&&(jcol>-1))
		{
			CHMAID[irow][jcol]=Stype; 
		}
		else
			throw std::logic_error("Error in SetCHMA");
			//int junk=1;

	}
	
	void SetCrop(int irow, int jcol, int ist, double per)
	{
		if ((irow<GetNumRows())&&(jcol<GetNumCols())&&(irow>-1)&&(jcol>-1) &&(ist<NStatus))
		{
			Trees[irow][jcol][ist]=per; 
			if ((per<0)||(per>1))
				throw std::logic_error("Error in SetCrop");
		}
		else
			throw std::logic_error("Error in SetCrop");
			//int junk=1;

	}

	
	void SetPlantaScore(int Plantacount,  double Score)
	{
		if (PlantationInfo[Plantacount].GetPlantNum() == Plantacount + 1)
			PlantationInfo[Plantacount].SetBlockScores(Score); // I think this is the plantation score
		else
			throw std::logic_error("PlantaCount and index don't match SPS2"); //the function should not exit here unless there is no Maize
	}


	void AddHealthLoss(int Plantacount, double incrInf)
	{
		if (PlantationInfo[Plantacount].GetPlantNum() == Plantacount + 1)
			PlantationInfo[Plantacount].AddHLoss(incrInf); // this adds to health loss should be reset every month
		else
			throw std::logic_error("PlantaCount and index don't match SPS1"); //the function should not exit here unless there is no Maize

	}

	void ResetHealthLoss(int Plantacount)
	{
		if (PlantationInfo[Plantacount].GetPlantNum() == Plantacount + 1)
			PlantationInfo[Plantacount].ResetHealthLoss(); // this adds to health loss should be reset every month
		else
			throw std::logic_error("PlantaCount and index don't match SPS3"); //the function should not exit here unless there is no Maize
		
	}

	double GetHealthLoss(int Plantacount)
	{
		if (PlantationInfo[Plantacount].GetPlantNum() == Plantacount + 1)
			return PlantationInfo[Plantacount].GetHLoss(); // this adds to health loss should be reset every month
		else
			throw std::logic_error("PlantaCount and index don't match SPS4"); //the function should not exit here unless there is no Maize

	}



	void SetPlantaControlBelief(int Plantacount, double myBelief)
	{
		if (PlantationInfo[Plantacount].GetPlantNum() == Plantacount + 1)
			PlantationInfo[Plantacount].SetPBelief(myBelief); // I think this is the plantation score
		else
			throw std::logic_error("PlantaCount and index don't match SPS5"); //the function should not exit here unless there is no Maize
	}

	void SetPlantaDisStatus(int PlantaID, int myStat)
	{
		int junk = PlantationInfo.size();
		if (PlantationInfo.size()>PlantaID)
		{
			if (PlantationInfo[PlantaID].GetPlantNum()==PlantaID+1)
				PlantationInfo[PlantaID].SetDisStatus(myStat);
			else
				throw std::logic_error("PlantaID and index don't match SPDS"); //the function should not exit here unless there is no Maize
		}
		else
			throw std::logic_error("PlantaID and index don't match SPDS2"); //the function should not exit here unless there is no Maize
	}
	

	int GetBlockCHMA(int BlockID)
	{
		return PlantationInfo[BlockID].GetFmCHMA();
	}

	int GetBlockPlanta(int BlockID)
	{
		return 0; //PlantationInfo[BlockID].GetGetPlanta();
	}

	double GetCrop(int irow, int jcol, int ist)
	{
		if ((irow>-1)&&(irow<Nrows)&&(jcol>-1)&&(jcol<Ncols))
			return Trees[irow][jcol][ist];
		else
			throw std::logic_error("error in get crop");
			//return -9;
	}

	int GetPlanta(int irow, int jcol)
	{
		if ((irow>-1)&&(irow<Nrows)&&(jcol>-1)&&(jcol<Ncols))
			return PlantaID[irow][jcol];
		else
			throw std::logic_error("error in get Planta");
			//return -9;
	}

	int GetCHMA(int irow, int jcol)
	{
		if ((irow>-1)&&(irow<Nrows)&&(jcol>-1)&&(jcol<Ncols))
			return CHMAID[irow][jcol];
		else
			throw std::logic_error("error in get CHMA");
			//int junk=-9;
	}

	int GetBlock(int irow, int jcol)
	{
		if ((irow>-1)&&(irow<Nrows)&&(jcol>-1)&&(jcol<Ncols))
			return BlockID[irow][jcol];
		else
			throw std::logic_error("error in get Block");
	}

	
	int GetMaxNumStatus(){return NStatus;}
			
	void GetOWLav(int irow, int jcol, int iyear, double& totHPop, double& totIPop){ GridBugs[irow][jcol].GetOWLav(iyear, totHPop, totIPop);}

	
	
	
	void ClearData()
	{
		for (int icount=0; icount<Nrows; icount++)
		{
			for (int jcount=0; jcount<Ncols; jcount++)
			{
				GridBugs[icount][jcount].Clear();
				for (int scount=0; scount<NStatus; scount++)
				{
					Trees[icount][jcount][scount]=0;
				}
			}
		}

		PlantationInfo.clear();

	}

	void ClearTempBugs()
	{
		for (int icount=0; icount<Nrows; icount++)
		{
			for (int jcount=0; jcount<Ncols; jcount++)
			{
				GridBugs[icount][jcount].ClearTempBug();
				
			}
		}
	}

	void SetTempBug(int irow, int jcol, int Health, double NewNum)
	{
		GridBugs[irow][jcol].SetTempBug(Health, NewNum);

	}
	double GetTempBug(int irow, int jcol, int Health)
	{
		if ((irow>-1)&&(irow<Grid::GetNumRows())&&(jcol>-1)&&(jcol<Grid::GetNumCols()))
			return GridBugs[irow][jcol].GetTempBug(Health);
		else
			throw std::logic_error("index out of bounds");
	}

	int GetControl(int irow, int jcol, int monthsAgo) //Get control 0 = no control 1 = control
	{
		if ((irow>-1)&&(irow<Grid::GetNumRows())&&(jcol>-1)&&(jcol<Grid::GetNumCols()))
			return G_Control[irow][jcol][monthsAgo];
		else
			throw std::logic_error("Error on Get Control index out of bounds");
	}

	void SetControl(int irow, int jcol, int Control) //Set control 0 = no control 1 = control
	{
		if ((irow>-1) && (irow<Grid::GetNumRows()) && (jcol>-1) && (jcol < Grid::GetNumCols()))
		{
			for (int icount = 0; icount < HistLen - 1; icount++)
			{
				G_Control[irow][jcol][icount + 1] = G_Control[irow][jcol][icount];
			}
			G_Control[irow][jcol][0] = Control;
		}
		
		else
			throw std::logic_error("Error in set control index out of bounds");
	}

	void SetIntVars(double jdist, double gridLen, double lambda)
	{
		IntVars[0]=jdist;
		IntVars[1]=gridLen;
		IntVars[2]=lambda;

	}

	void GetIntVars(double& jdist, double& gridLen, double& lambda)
	{
		jdist=IntVars[0];
		gridLen=IntVars[1];
		lambda=IntVars[2];
		
	}

	
	void SetInvasion(int irow, int jcol, int yr) //stores where Bugs have reached in an invasion from left to right
	{
		invasionYr[irow][jcol]=yr;
	}

	int GetInvasion(int irow, int jcol)
	{
		return invasionYr[irow][jcol];
	}

	void SetLand(int newL, double alpha_W, double beta_W, double alpha_B, double beta_B, int& iflag) //double alpha_W, double beta_W, double alpha_B, double beta_B are paraters for beta dristribution describing initial belief structures 
	{
		
		int UBlockID=1; //unique Block ID	

		for (int icount=0;icount<Nrows; icount++)
		{
			for (int jcount=0;jcount<Ncols; jcount++)
			{
				SetCHMA(icount,jcount,-9); // 0 =minnesota and 1= wisconsin
				SetPlanta(icount,jcount,-9);
				SetBlock(icount,jcount,-9);
				for (int kcount=0; kcount<NStatus; kcount++)
					SetCrop(icount,jcount,kcount,0);  //set all status to zero
				
				SetCrop(icount,jcount,NStatus-1,1);  //set non-trees to 100%

			}
		}
		if (newL==0) //then read grid spec from file (not coded yet)
		{
			char path_buffer[_MAX_PATH];
			char Npath_buffer[_MAX_PATH];
			_getcwd( path_buffer, _MAX_PATH );
			strcpy_s(Npath_buffer, path_buffer);
		//	strcat_s(Npath_buffer, "\\CHMA15\\Plant15.txt"); //was counties
		//	strcat_s(Npath_buffer, "\\CHMA54\\Plant54.txt"); //was counties
			strcat_s(Npath_buffer, "\\CHMAIR\\PlantIR.txt"); //was counties grid size rows=350 cols = 400
			std::ifstream    InF(Npath_buffer);

			strcpy_s(Npath_buffer, path_buffer);
		//	strcat_s(Npath_buffer, "\\CHMA15\\Block15.txt"); //was Blocks
		//	strcat_s(Npath_buffer, "\\CHMA54\\Block54.txt"); //was Blocks
			strcat_s(Npath_buffer, "\\CHMAIR\\PlantIR.txt"); //was counties grid size rows=350 cols = 400
			std::ifstream    InFF(Npath_buffer);

			strcpy_s(Npath_buffer, path_buffer);
		//	strcat_s(Npath_buffer, "\\CHMA15\\CHMA15.txt"); //was CHMAs
		//	strcat_s(Npath_buffer, "\\CHMA54\\CHMA54.txt"); //was Blocks
			strcat_s(Npath_buffer, "\\CHMAIR\\PlantIR.txt"); //was counties grid size rows=350 cols = 400
			std::ifstream	 InFFF(Npath_buffer);

			/*strcpy_s(Npath_buffer, path_buffer);
			strcat_s(Npath_buffer, "\\CHMA15\\TreesIni15.txt"); //was crops
			std::ifstream	 InC(Npath_buffer);*/

			int TempPlanta, TempBlock, TempCHMAID;
			int Sick(0); //Flag to see of sick tree set/
			int iyr=2;
			int MaxPlantID=-1; // This keeps track of the maximum ID number for the plantations for use later 

			if (InF&&InFF&&InFFF)
			{
				
				for (int jcount=0; jcount<Nrows; jcount++)
				{
					for (int icount=0; icount<Ncols; icount++)
					{
						InF>>TempPlanta; 
						PlantaID[jcount][icount]=TempPlanta; 
						if (TempPlanta>MaxPlantID)
							MaxPlantID=TempPlanta;
						InFF>>TempBlock;
						BlockID[jcount][icount]=TempBlock;
						InFFF>>TempCHMAID;
						CHMAID[jcount][icount]=TempCHMAID;
												
						//Here we set the first cell that has a crop to 10% inf not syp
						if (TempBlock!=0)
						{
							if (Sick<3)
							{
								SetCrop(jcount, icount, 0, 0.9); //All healthy
								SetCrop(jcount, icount, 1, 0.1); //No Infected not syp
								SetCrop(jcount, icount, 2, 0.0); //No Infected 
								SetCrop(jcount, icount, 3, 0); //No symp
								SetCrop(jcount, icount, 4, 0); //No removed
								Sick++;
							}
							else
							{
								SetCrop(jcount, icount, 0, 1.0); //All healthy
								SetCrop(jcount, icount, 1, 0); //No Infected not syp
								SetCrop(jcount, icount, 2, 0); //No Infected 
								SetCrop(jcount, icount, 3, 0); //No symp
								SetCrop(jcount, icount, 4, 0); //No removed
								
							}

						}
						else
						{
							SetCrop(jcount, icount, 0, 0); //No healthy
							SetCrop(jcount, icount, 1, 0); //No Infected not syp
							SetCrop(jcount, icount, 2, 0); //No Infected 
							SetCrop(jcount, icount, 3, 0); //No symp
							SetCrop(jcount, icount, 4, 1); //All  removed

						}
							
					}
					
				}
				InF.close();
				InFF.close();
				InFFF.close();
				//InC.close();
				

			}
			else
				throw std::logic_error("Error opening inpur files");

			Grid::SetTotalPlantasFilled(-1);
			int mySymp(0); //Set to zero is max syptomatic trees in any one cell in the plantaion is less than 0.05
			double myWorry(0); //start unworried
			 
			int myCHMA(0);
			
			for (int pcount=0; pcount<MaxPlantID; pcount++)
			{
				//int myCells[MaxNumPCells]; 
				std::vector<int> myCells;
				int countcell=0;
				for (int jcount=0; jcount<Nrows; jcount++)
				{
					for (int icount=0; icount<Ncols; icount++)
					{
						if (pcount+1==PlantaID[jcount][icount])
						{
							myCells.push_back(CoordToGridID(jcount, icount)); //gathering the co-ordinates of all cells in a plantation(grove)
							myCHMA=CHMAID[jcount][icount];
							
						}
					}
				}
				if (myCells.size()<1)
					throw std::logic_error("Number of cells in a plantation is less that 1");

				//generate some random belief in control strategy//////////////
				int ifail = 1;
				int lstate(0);
				int subid(1);
				int *state;
				state = new int[lstate];
				int Lseed = 1;
				int seed[1];
				seed[0] = 1;
				G05KGF(1, subid, state, lstate, ifail);
				//G05KFF(1,subid, seed, Lseed, state,lstate,ifail); //always get same vale of random permutaion for a given seed
				delete[] state;
				state = new int[lstate];
				ifail = 1;
				G05KGF(1, subid, state, lstate, ifail);
				//G05KFF(1,subid, seed, Lseed, state,lstate,ifail);
				double ubeta[1]; //first is probability that we set a bug, second and third define which cell in which plant we choose
				ifail = 1;
		//		G05SAF(1, state, urand, ifail);
				G05SBF(1, alpha_W, beta_W, state, ubeta, ifail);
				myWorry = ubeta[0]; //SET INITAIL WORRY

				double myBelief(0);

				//If the random number is larger than propbelief we don't believe
				G05SBF(1, alpha_B, beta_B, state, ubeta, ifail);
				myBelief = ubeta[0]; //SET INITAIL WORRY
				/////////////////////////////////////////////////////
			

				PInfo TempFInfo(myCHMA, pcount+1, mySymp, myWorry, myBelief, myCells);
				
				
				PlantationInfo.push_back(TempFInfo);
			}

		}
		else //Set up a defind simple structure
		{
			int btrip=1; //keeps track of each block length
			int blockN=-1; //block number will start at zero
			for (int icol=0; icol<Ncols; icol++)
			{
				for (int jrow=0; jrow<Nrows; jrow++)
				{
					
					//allocate blocks of size 4x1
					if ((btrip>4)||(jrow==0)) //block number jumps after 4 or at end of column
					{
						blockN=blockN+1;
						btrip=1; //keeps track of each block length
					}
					Grid::SetBlock(jrow, icol, blockN);
					for (int ist=1;ist< NStatus; ist++)
						Grid::SetCrop(jrow,icol,ist, 0); // Set infected, sympomatic and removed to zero
					Grid::SetCrop(jrow,icol,0, 1); // Set healthy to 1
					btrip++;
				}
			} 


			//Fit in CHMA there are 2x4 plantations in a CHMA

			int Rstart=0;
			int Cstart=0;
			int CColN=30*4; //There are 6x30 blocks in a plantation and our blocks are 4x1 in size. There are 2x4 plantations in a block This is the space they will take
			int CRowN=6*4*2; //There are 6x30 blocks in a plantation and our blocks are 4x1 in size. 
			int DoneFlag=0;
			int CHMAN=-1; //block number will start at zero
			if(CColN>Ncols)
					CColN=Ncols;
			if(CRowN>Nrows)
					CRowN=Nrows;
			while (DoneFlag==0)
			{
				CHMAN=CHMAN+1;
				for (int icol=Cstart; icol<CColN; icol++)
				{
					for (int jrow=Rstart; jrow<CRowN; jrow++)
					{	
						Grid::SetCHMA(jrow, icol, CHMAN);
					}
				} 
				//Keep moving across column blocks
				Cstart=CColN;
				CColN=CColN+30*4;
				if(CColN>Ncols)
					CColN=Ncols;
				if (Cstart==Ncols) //We need to go down a chunk
				{
					Cstart=0;
					CColN=30*4;
					if (CColN>Ncols)
						CColN=Ncols;
					
					Rstart=CRowN;
					CRowN=CRowN+6*4*2;
					if (CRowN>Nrows)
						CRowN=Nrows;
					if (Rstart==Nrows) //We need to stop
					{
						DoneFlag=1;
					}
				}


			}

			//Fit plantations in last so we can allocate plantation struvture
			Rstart=0;
			Cstart=0;
			int PColN=30; //There are 6x30 blocks in a plantation and our blocks are 4x1 in size. This is the space they will take
			int PRowN=6*4; //There are 6x30 blocks in a plantation and our blocks are 4x1 in size. 
			DoneFlag=0;
			int PlantN=0; //block number will start at one
			if(PColN>Ncols)
					PColN=Ncols;
			if(PRowN>Nrows)
					PRowN=Nrows;
			while (DoneFlag==0)
			{
				PlantN=PlantN+1;
				int myCHMA(0);
				
				//int myCells[MaxNumPCells]; 
				std::vector<int> myCells;
				for (int icol=Cstart; icol<PColN; icol++)
				{
					for (int jrow=Rstart; jrow<PRowN; jrow++)
					{	
						Grid::SetPlanta(jrow, icol, PlantN);
						myCHMA=Grid::GetCHMA(jrow, icol);
						myCells.push_back(CoordToGridID(jrow, icol));
						

					}
				} 

				
				int mySymp(0); //Set to zero is max syptomatic trees in any one cell in the plantaion is less than 0.05
				double myWorry(0); //worry about possible infection 
				double myBelief(0); 

				//generate some random belief in control strategy
				int ifail = 1;
				int lstate(0);
				int subid(1);
				int *state;
				state = new int[lstate];
				int Lseed = 1;
				int seed[1];
				seed[0] = 1;
				G05KGF(1, subid, state, lstate, ifail);
				//G05KFF(1,subid, seed, Lseed, state,lstate,ifail); //always get same vale of random permutaion for a given seed
				delete[] state;
				state = new int[lstate];
				ifail = 1;
				G05KGF(1, subid, state, lstate, ifail);
				//G05KFF(1,subid, seed, Lseed, state,lstate,ifail);
				double ubeta[1]; //first is probability that we set a bug, second and third define which cell in which plant we choose
				ifail = 1;
				//		G05SAF(1, state, urand, ifail);
				G05SBF(1, alpha_W, beta_W, state, ubeta, ifail);
				myWorry = ubeta[0]; //SET INITAIL WORRY

				myBelief=0;

				//If the random number is larger than propbelief we don't believe
				G05SBF(1, alpha_B, beta_B, state, ubeta, ifail);
				myBelief = ubeta[0]; //SET INITAIL WORRY
				/////////////////////////////////////////////////////



				PInfo TempFInfo(myCHMA, PlantN, mySymp, myWorry, myBelief, myCells);
				PlantationInfo.push_back(TempFInfo);
				if (PlantationInfo[PlantN-1].GetPlantNum()!=PlantN)			
					throw std::logic_error("Error in Farm indexing"); //the function should not exit here unless there is no Maize
				
			
				

				//Keep moving across column blocks
				Cstart=PColN;
				PColN=PColN+30;
				if(PColN>Ncols)
					PColN=Ncols;
				if (Cstart==Ncols) //We need to go down a chunk
				{
					Cstart=0;
					PColN=30;
					if (PColN>Ncols)
						PColN=Ncols;
					Rstart=PRowN;
					PRowN=PRowN+6*4;
					if (PRowN>Nrows)
						PRowN=Nrows;
					if (Rstart==Nrows) //We need to stop
					{
						DoneFlag=1;
					}
				}


			}

			


			}


		char path_buffer[_MAX_PATH];
		char Npath_buffer[_MAX_PATH];
		_getcwd( path_buffer, _MAX_PATH );
		strcpy_s(Npath_buffer, path_buffer);

		strcpy_s(Npath_buffer, path_buffer);
	    strcat_s(Npath_buffer, "\\OutFiles\\Plantas.txt");
		std::ofstream    OutF;
		OutF.open(Npath_buffer, std::ios_base::out);

		if (OutF)
			int junk=1;


		strcpy_s(Npath_buffer, path_buffer);
	    strcat_s(Npath_buffer, "\\OutFiles\\Blocks.txt");
		std::ofstream    OutFF;
		OutFF.open(Npath_buffer, std::ios_base::out);

		strcpy_s(Npath_buffer, path_buffer);
	    strcat_s(Npath_buffer, "\\OutFiles\\CHMAs.txt");
		std::ofstream	 OutFFF;
		OutFFF.open(Npath_buffer, std::ios_base::out);

		

		
		for (int jcount=0; jcount<Nrows; jcount++)
		{
			for (int icount=0; icount<Ncols; icount++)
			{
				OutF<<GetPlanta(jcount,icount)<<' ';
				OutFF<<GetBlock(jcount,icount)<<' ';
				OutFFF<<GetCHMA(jcount,icount)<<' ';
				
				OutF.flush();
				//OutC.flush();
				
			}
			OutF<<'\n';
			OutF.flush();
			OutFF<<'\n';
			OutFFF<<'\n';
			//OutC<<'\n';
			//OutC.flush();
			
		}
		/*if (OutC.bad()) {
			int i = 0;
		}*/
		OutF.close();
		OutFF.close();
		OutFFF.close();
		//OutC.close();
	

	} //end of function SetLand

	


}
////////////////////////////END OF NAMESPACE GRID///////////////////////////////////////////

namespace FlightP // keeps a vector of the proportions of Bugs of each type that fly each distance
	
{
	std::vector<double> BugFlight; //
	
	void ClearData()
	{
		
		BugFlight.clear(); //
		

	};


	void MakeProportions(std::vector<double>& Type, double lambda)
	{
			//Type 1 (Type[0]) centre cell for other types see 
			int idist=0;
			int jdist=0;
			
			Type.push_back(Integrate(idist, jdist, lambda)); //position 0
			Type.push_back(Integrate(0, 1,  lambda)); //position 1
			Type.push_back(Integrate(1, 1, lambda)); //position 2
			//Test to see how much of the population is in this inner square
			double popT=Type[0]+4*Type[1]+4*Type[2];
			//if (popT<0.99) //keep calculating
			int Mcount=3;
			while (popT<0.975) 
			{
				for (int icount=0; icount<Mcount; icount++)
				{
					double Ival=Integrate(icount, Mcount-1, lambda);
					Type.push_back(Ival);
					if ((icount>0)&&(icount<Mcount-1))
						popT=popT+8*Ival;
					else
						popT=popT+4*Ival;
				}
				Mcount++;
				if (Mcount>900)
					int junk=0;
								
			}
			//put remained in the centre
			
			Type[0]=Type[0]+1.0-popT;
			if (Type[0]<0) Type[0]=0;
	

	}

	void Initialise()
	{ 
			double lambda = 0.00035; // For weekly time step. See page 41 of notebook
			MakeProportions(BugFlight, lambda);
		


	}

	void GetBugFlight(std::vector<double>& MyDist)
	{
		int Len=BugFlight.size();
		for (int icount=0; icount<Len; icount++)
		{
			MyDist.push_back(BugFlight[icount]);
		}
	}

	

		
} //end of namespace

double Integrate(int idist, int jdist, double lambda)
{
	//This function integrates the exponential distributions over each cell to see how many bugs we have in a cell.

	//idist= irow-i0 where irow is place we are estimating population and i0 is where Buge started
	//jdist=jrow-j0 ditto
	double gridLen=100.0; // grid length is 100 m
	
	Grid::SetIntVars(jdist, gridLen, lambda);

	double xlb[2], xub[2];
	xlb[0]=(idist -0.5)*gridLen;
	xlb[1]=(jdist -0.5)*gridLen;
	xub[0]=(idist +0.5)*gridLen;
	xub[1]=(jdist +0.5)*gridLen;

	

	
	double ABSACC=0.0;
	double RELACC=0.0000001;
	//double ANS(0);
	int NDIM(2);
	int NumFunc(1);
	int MaxCals(180000);
	int MinCals(10);
	const int LENWRK=50000;//(NDIM+NumFunc+2)*(10+MaxCals);
	double WRKSTR[LENWRK];
	int IFAIL=1;
	double ANS[1];
	double ABSET[1];
		
	
	//D01DAF(xlb, xub, ylbFunc, yubFunc,fFunc,ABSACC, ANS, NPTS, IFAIL);
	D01EAF(NDIM, xlb, xub, MinCals, MaxCals, NumFunc, fFunc2, ABSACC, RELACC, LENWRK, WRKSTR, ANS, ABSET, IFAIL);

	if (ANS[0]<0)
		throw std::logic_error("Error in integration");
	if (IFAIL!=0)
		int junk=0;
	
	//ANS[0]=ANS[0]*ga;
	return ANS[0];

}


extern "C" void __stdcall fFunc2(const int& NDIM, const double x[], const int& NumFunc, double funcY[] )
//void  fFunc2(const int& NDIM, const double x[], const int& NumFunc, double funcY[] )
{
	double jdist, gridLen, lambda;
	Grid::GetIntVars(jdist, gridLen, lambda);
	int IFAIL=1;
	double r=pow(x[0],2.0)+pow(x[1],2.0);
	r=pow(r, 0.5);
	double func(0);
	double OOtwoPi=1.0/(2.0*3.141572);
	if (r>0)
	{
		//func=OOtwoPi*lambda*lambda*exp(-lambda*r);
		func=OOtwoPi*lambda*exp(-lambda*r)/r;
	}
	else
	{
		func=0;
	}

	funcY[0]=func;


}

void GridIDToCoord(int Floc, int& irow, int& icol)
{
	irow = Floc % Grid::GetNumRows();  // % is remainder after dividing
	icol=(floor(Floc/Grid::GetNumRows())); 
}
int CoordToGridID(int jrow, int icol)
{
	return icol*Grid::GetNumRows()+jrow;
}

