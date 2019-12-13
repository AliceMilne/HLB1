#include "HLBGrid.h"
#include "HLBModel.h"
#include <math.h>
#include <fstream>
#include <iostream>
#include <direct.h>
#include <list>
//#include "d:\\Program Files (x86)\\NAG\\FL25\\fldll254ml\\c_headers\\nagmk25.h"

// Parameters
double k_cap = 150000; //carrying capacity. 1500 trees per block 300 bugs per tree - needs revising 3 cells in a block 
double sigma= 12; //number of off spring per vector at low density - check with Frank
double infRate =  0.0033; //  infection rate of bug ... 
double delta_inf=1.0; //between 0 and 1 how infectius 'infected' is compared with 'symptomatic'

double Talpha = 0.0033; // infRate; //now assuming same as alpha but was 0.01; //probability a tree gets infected given Y infected vectors. Alaph*Y = mean number of infections
double Trho=1; // 1/gamma = mean time that trees are in infected but not infectious state Now two months
double Tgamma=0.1666; // 1/gamma = mean time to go from state I to state S

//0=nontransgenic maize; 1=transgenic wheat maize; 2=RIB; 3=non-prefferential host; 4=uninhabitable area; 
double const RRCropSurv[5]={1, 1, 1, 1, 0}; //survival adjustment due to cropping
double const SSCropSurv[5]={1, 0.001, 0.001, 1, 0.0}; //survival adjustment due to cropping

double WorryThreshold=0.002; // this equates to a 1 tree a cell


void InitialiseModel(int method, double alpha_W, double beta_W, double alpha_B, double beta_B)
{
	
	Grid::ClearData();
	FlightP::ClearData();
//Initialises some outputfiles that we shall write our results in. 
	//Currently commented out as we don't know what we want to write out
	
	char path_buffer[_MAX_PATH];
	char Npath_buffer[_MAX_PATH];
	_getcwd( path_buffer, _MAX_PATH );
	strcpy_s(Npath_buffer, path_buffer);
	strcat_s(Npath_buffer, "\\OutFiles");
	int  dirMade=_mkdir(Npath_buffer);

/*	strcpy_s(Npath_buffer, path_buffer);
	strcat_s(Npath_buffer, "\\OutFiles\\Population1.txt");
	std::ofstream File1(Npath_buffer);

	strcpy_s(Npath_buffer, path_buffer);
	strcat_s(Npath_buffer, "\\OutFiles\\Population2.txt");
	std::ofstream File2(Npath_buffer);

	strcpy_s(Npath_buffer, path_buffer);
	strcat_s(Npath_buffer, "\\OutFiles\\Population3.txt");
	std::ofstream File3(Npath_buffer);

	strcpy_s(Npath_buffer, path_buffer);
	strcat_s(Npath_buffer, "\\OutFiles\\Population4.txt");
	std::ofstream File4(Npath_buffer);

	strcpy_s(Npath_buffer, path_buffer);
	strcat_s(Npath_buffer, "\\OutFiles\\Population5.txt");
	std::ofstream File5(Npath_buffer);

	strcpy_s(Npath_buffer, path_buffer);
	strcat_s(Npath_buffer, "\\OutFiles\\Population6.txt");
	std::ofstream File6(Npath_buffer);
	
	File1.close();
	File2.close();
	File3.close();
	File4.close();
	File5.close();
	File6.close();


	*/
	
	if (method==0) //	Sets simple croping of all trees
	{
		for (int icount=0; icount<Grid::GetNumRows(); icount++)
		{
			for (int jcount=0; jcount<Grid::GetNumCols(); jcount++)
			{
			    int myStatus=0; //Status of trees. This is either Healthy=0, 1=Infected not infectious, 2=Infectious no symtoms=1, Symptoms=2, Removed (which includes the area never planted) =4 
				for (int ist=1;ist<Grid::GetMaxNumStatus(); ist++)
						Grid::SetCrop(icount,jcount,ist, 0); // Set infected, sympomatic and removed to zero
				Grid::SetCrop(icount, jcount, myStatus, 1);
				Grid::SetInvasion(icount, jcount, -1);
				Grid::SetControl(icount, jcount, 0);

				
			}
		}

		
		
		Grid::SetCrop(100, 100, 0, 0.988); // Set healthy to 0.988
		Grid::SetCrop(100, 100, 3, 0.002); // Set sympo to 0.002
		Grid::SetCrop(100, 100, 2, 0.005); // Set inf   to 0.005
		Grid::SetCrop(100, 100, 1, 0.005); // Set inf not infrctious to 0.5
		

	}
	else if(method==1) //reads from landscape files. This is used in the paper
	{

		SetLand(0, alpha_W, beta_W, alpha_B, beta_B); //0 = read from files
			for (int icount=0; icount<Grid::GetNumRows(); icount++)
			{
				for (int jcount=0; jcount<Grid::GetNumCols(); jcount++)
				{
					Grid::SetInvasion(icount, jcount, -1);
				}
			}
	}
	else //sets up some blocking stucture but continual trees
	{
		SetLand(1, alpha_W, beta_W, alpha_B, beta_B);
			for (int icount=0; icount<Grid::GetNumRows(); icount++)
			{
				for (int jcount=0; jcount<Grid::GetNumCols(); jcount++)
				{
					Grid::SetInvasion(icount, jcount, -1);
				}
			}

			
			Grid::SetCrop(100, 100, 0, 0.988); // Set healthy to 0.5
			Grid::SetCrop(100, 100, 3, 0.002); // Set sympo to 0.5
			Grid::SetCrop(100, 100, 2, 0.005); // Set inf   to 0.5
			Grid::SetCrop(100, 100, 1, 0.005); // Set inf not infrctious to 0.5


	}

	
	for (int icount=0; icount<Grid::GetNumRows(); icount++)
	{
		for (int jcount=0; jcount<Grid::GetNumCols(); jcount++)
		{
			if (Grid::GetCrop(icount, jcount,4) < 1) ////Status of trees. Status 4 is 'removed which includes never planted there. If this is the case there will be no bug
			{

				
				Grid::SetNumG(icount, jcount, 0, 50000.0); //Sets number healthy bugs
				Grid::SetNumG(icount, jcount, 1, 0.0); //Sets number of infected bugs
			
				
			}
			else
			{
				Grid::SetNumG(icount, jcount, 0, 0); //Sets number healthy bugs
				Grid::SetNumG(icount, jcount, 1, 0); //Sets number of infected bugs
				
			}
			
		
		}
	}

	

	
	FlightP::Initialise();
	
	



}

void DeleteModel()
{
	Grid::ClearData();
	FlightP::ClearData();
	

}

void SetLand(int readin, double alpha_W, double beta_W, double alpha_B, double beta_B)
{
	int iflag(0);
	Grid::SetLand(readin, alpha_W, beta_W, alpha_B,  beta_B, iflag);
	if (iflag==1)
		throw std::logic_error("Error in Planta assignment");
	else if (iflag==2)
		throw std::logic_error("Error creating Planta landscape");
	else if (iflag==3)
		throw std::logic_error("Error creating Block assignment");
						

}

//Eggs are laid on tips of growing shoots on and between unfurling leaves. Females lay 300 to 800 eggs during their lifetime. 
//Nymphs pass through five instars. The total life cycle requires from 15 to 47 days, depending on environmental factors such as temperature and season. 
//The adults may live for more than a month. There is no diapause, but populations are typically low in the winter or during dry periods. 
//There are 9 to 10 generations a year, with up to 16 observed under observation in field cages. 


void ModelMonth(int imth, double randInf, double kill_rate, int frequControl)
{

	/*char path_buffer[_MAX_PATH];
	char Npath_buffer[_MAX_PATH];
	_getcwd( path_buffer, _MAX_PATH );
	strcpy_s(Npath_buffer, path_buffer);
	strcat_s(Npath_buffer, "\\OutFiles\\MyErrors.txt");
	
	std::ofstream ofm(Npath_buffer, std::ofstream::out | std::ofstream::app);

	ofm<<imth<<'\n';
	ofm.flush();*/

	for (int icount=0; icount<Grid::GetNumRows(); icount++)
	{
		for (int jcount=0; jcount<Grid::GetNumCols(); jcount++)
		{
			
			Generation(icount, jcount, imth, kill_rate,  frequControl) ; //grows bugs and calculates the number of ingected ones based on tree health
		}
	}

//	WriteBugs(imth); // this function writes a map of infected bugs for total see Tbugs

	Flight(); //Bugs distribute
	Flight(); //Bugs distribute
	Flight(); //Bugs distribute
	Flight(); //Bugs distribute
	
	if (randInf>0)
	{

		if (randInf>1)
		{
			throw std::logic_error("probability must be less than 1");
		}
		else
		{

			int myNumPlantas=Grid::GetPlantaInfoSize();
			//pick one at random
			int ifail=1;
			int lstate(0);
			int subid(1);
			int *state;
			state=new int[lstate];
			int Lseed=1;
			int seed[1];
			seed[0]=1;
			G05KGF(1,subid,state,lstate,ifail);
			//G05KFF(1,subid, seed, Lseed, state,lstate,ifail); //always get same vale of random permutaion for a given seed
			delete [] state ;
			state=new int[lstate];
			ifail=1;
			G05KGF(1,subid,state,lstate,ifail);
			//G05KFF(1,subid, seed, Lseed, state,lstate,ifail);
			double urand[3]; //first is probability that we set a bug, second and third define which cell in which plant we choose
			ifail=1;
			G05SAF (3,  state, urand, ifail);

			if(urand[0]<randInf)
			{
				int myPlanta=int(floor(urand[1]*myNumPlantas));
			    //how many cells in this planta
				int myNumCells=Grid::GetNumCells(myPlanta);
				//Randomly select one
				int Ccount=floor(urand[2]*myNumCells);
				int Floc=Grid::GetCellLocation(myPlanta, Ccount); //this returns the code for the index of the cell
				int irow, icol;
				GridIDToCoord(Floc, irow, icol); //converts grid id to gris coordinates
				//Set infected bug =1 if not already some there.
				Grid::SetNumG(irow, icol, 1, 1);
			}

		}

	}



	for (int icount=0; icount<Grid::GetNumRows(); icount++)
	{
		for (int jcount=0; jcount<Grid::GetNumCols(); jcount++)
		{
			if ((icount == 10) && (jcount == 343))
				int junk = 0;
			Tree_Generation(icount, jcount, imth); //Bugs land on trees change health status and things proceed
			
		}
	}

	

			
}

void Generation(int irow, int jcol, int imth, double kill_rate, int frequControl) // Model the population dynamics over a month 
{
	if ((irow == 1)&(jcol == 152))
		int junk = 0;

	char path_buffer[_MAX_PATH];
	char Npath_buffer[_MAX_PATH];
	_getcwd( path_buffer, _MAX_PATH );
	strcpy_s(Npath_buffer, path_buffer);
	strcat_s(Npath_buffer, "\\OutFiles\\MyErrors.txt");
	
//	std::ofstream ofm(Npath_buffer, std::ofstream::out | std::ofstream::app);



	double Bug_H=0;
	double Bug_I=0;
	
	//Get the eggs laid by healthy females 
	GetPop(irow, jcol, Bug_H, Bug_I);  
	if (Bug_H+Bug_I>0)
	{
		
			//Increase in bugs per time step
			double spray = 0.0; // If they buy into control strategy and its the month when control is applied then spray
			if (Grid::GetControl(irow, jcol, 0) > 0)
			{
				double denom = double(frequControl) / 12.0;
				double ans1 = double(imth)*denom;
				ans1 = ans1 - floor(ans1);
				
				if (ans1<0.0000001)
					spray = 1.0;
			}
			


			double theta = spray*kill_rate; // control applied to current time step
			double TBug = Bug_H + Bug_I;
			TBug = k_cap*TBug / (k_cap / sigma + TBug)*(1 - theta);
			//Bug_I=k_cap*Bug_H/(k_cap/sigma+TBug)*(1-theta);
			double propInf = Grid::GetCrop(irow, jcol, 2); //proprtion of infected trees and infectious. //Status of trees. This is either Healthy=0, Infected no symtoms=1, Symptoms=2, Removed (which includes the area never planted) =3 
			double propSymp = Grid::GetCrop(irow, jcol, 3); //proprtion of symptomatic

			// a proportion of bugs get infected
			Bug_I = (1 - exp(-infRate*(delta_inf*propInf + propSymp)))*TBug;

			Bug_I = round(Bug_I);
			Bug_H = TBug - Bug_I;

			Bug_H = round(Bug_H); 
			if (Bug_H<0) 
				Bug_H=0.0;

			Grid::SetNumG(irow, jcol, 0, Bug_H);
			Grid::SetNumG(irow, jcol, 1, Bug_I);
		
		
		

	

	}
	else
	{
		double Bug_T=0;
		Grid::SetNumG(irow, jcol, 0, Bug_T);
		Grid::SetNumG(irow, jcol, 1, Bug_T);
	}

	//ofm.close();
}



void Tree_Generation(int irow, int jcol, int iyr) //progesses tree health status in a time step
{

	char path_buffer[_MAX_PATH];
	char Npath_buffer[_MAX_PATH];
	_getcwd( path_buffer, _MAX_PATH );
	strcpy_s(Npath_buffer, path_buffer);
	strcat_s(Npath_buffer, "\\OutFiles\\MyErrors.txt");
	
		
	double Bug_H=0;
	double Bug_I=0;
	
	//Get the eggs laid by healthy females that are RR etc in this grid cell
	GetPop(irow, jcol, Bug_H, Bug_I);  
	if (Bug_I > 0)
		int junk = 0;
	double H_trees=Grid::GetCrop(irow, jcol, 0); //Healthy trees
	double E_trees=Grid::GetCrop(irow, jcol, 1); // Infected not infectious
	double I_trees=Grid::GetCrop(irow, jcol, 2); // Infected not symptomatic
	double S_trees=Grid::GetCrop(irow, jcol, 3); // Symptomatic



	
	double NH_trees=H_trees*exp(-Talpha*Bug_I);
	double NE_trees=E_trees*(1-Tgamma)+H_trees*(1-exp(-Talpha*Bug_I));
	double NI_trees=Tgamma*E_trees+I_trees*(1-Trho); 
	double NS_trees=Trho*I_trees+S_trees; 
	
	//If the number of susceptible trees has increased then count the added loss of health.

	double incrInf = NS_trees - S_trees;
	int myPlanta = Grid::GetPlanta(irow, jcol);
	if (incrInf>0)
	{
		Grid::AddHealthLoss(myPlanta-1, incrInf);
	}


	if (NS_trees>WorryThreshold) //if the number of syptimatic trees is > than a certain level
	{
		Grid::SetPlantaDisStatus(myPlanta-1, 1); //planta number -1 is the index in the vector

	}


	Grid::SetCrop(irow, jcol, 0, NH_trees);
	Grid::SetCrop(irow, jcol, 1, NE_trees);
	Grid::SetCrop(irow, jcol, 2, NI_trees);
	Grid::SetCrop(irow, jcol, 3, NS_trees);
	


}

void GetPop(int irow, int jcol, double& Bug_H, double& Bug_I )
{
	Bug_H=Grid::GetNumG(irow, jcol, 0); //get healthy bugs
	Bug_I=Grid::GetNumG(irow, jcol, 1); //get infected bugs
	
}



void Flight()
{
	//Disperse first
	//Clear Temp Vector
	Grid::ClearTempBugs();
	int season(1); //summer
	int post(0); //pre=0 post=1
	
//	WriteBugsT(1);
	
	//Bugs Fly 
	for (int icount=0; icount<Grid::GetNumRows(); icount++)
	{
		for (int jcount=0; jcount<Grid::GetNumCols(); jcount++)
		{
				double H_Bug=Grid::GetNumG(icount, jcount ,0); //healthy 
			//	Distribute2(icount, jcount, 0, H_Bug); //ALICE SWITCH ON BELOW BIT TOO LABELLED ALICE
				
				double I_Bug=Grid::GetNumG(icount, jcount, 1); //infected
				if (I_Bug>0)
					Distribute2(icount, jcount, 1, I_Bug);
		}
	}
	//Set grid of Bugs with new values 


	for (int icount=0; icount<Grid::GetNumRows(); icount++)
	{
		for (int jcount=0; jcount<Grid::GetNumCols(); jcount++)
		{
			//	Grid::SetNumG(icount, jcount, 0, Grid::GetTempBug(icount, jcount, 0)); //healthy  ALICE SWITCH BACK ON IF FLY HEALTHY
				Grid::SetNumG(icount, jcount, 1, Grid::GetTempBug(icount, jcount, 1));  //infected
			//	Grid::SetOWLav(icount, jcount, Grid::GetTempBug(icount, jcount, 0), Grid::GetTempBug(icount, jcount, 1));
				Grid::SetOWLav(icount, jcount, Grid::GetNumG(icount, jcount, 0), Grid::GetTempBug(icount, jcount, 1));
		}
	}
	
	
	//WriteBugsT(2);
	
}

void SetControl(int imth, double Per) //sets control across landscape in month imnth. If Per is negative then control is read in from files. If positive the proprtion is randomly distibuted across cells 
{
	//clear control
	for (int icount=0; icount<Grid::GetNumRows(); icount++)
	{
		for (int jcount=0; jcount<Grid::GetNumCols(); jcount++)
		{
				Grid::SetControl(icount, jcount, 0);
		}
	}
	
	if (Per<0)
	{
		//TODO read from files
	}
	else if (Per<=1)
	{
		
		//There are N plantations numbered 0 to n-1. 
		int Nplanta=Grid::GetPlantaInfoSize();
		int *INDX;
		INDX=new int[Nplanta];
		for (int icount=0; icount<Nplanta; icount++)
		{
			INDX[icount]=icount; // this is the vector we shall randomly permutate
		}

		int ifail=1;
		int lstate(0);
		int subid(1);
		int *state;
		state=new int[lstate];
		int Lseed=1;
		int seed[1];
		seed[0]=1;
		G05KGF(1,subid,state,lstate,ifail);
		//G05KFF(1,subid, seed, Lseed, state,lstate,ifail); //always get same vale of random permutaion for a given seed
		delete [] state ;
		state=new int[lstate];
		ifail=1;
		G05KGF(1,subid,state,lstate,ifail);
		//G05KFF(1,subid, seed, Lseed, state,lstate,ifail);
		double mu=11.9455;
		double sigsq= 0.2810*0.2810;

		ifail=1;
		G05NCF (INDX, Nplanta, state, ifail);

		//Now the first Per percent will spray
		int nspray=floor(Per*Nplanta); // the number that spray
		for (int icount=0; icount<nspray; icount++)
		{
			int myPlanta=INDX[icount];
			int myNumCells=Grid::GetNumCells(myPlanta);
			for (int jcount=0; jcount<myNumCells; jcount++)
			{
				int Floc=Grid::GetCellLocation(myPlanta, jcount); //this returns the code for the index of the cell
				int irow, icol;
				GridIDToCoord(Floc, irow, icol); //converts grid id to gris coordinates
				Grid::SetControl(irow, icol, 1);
			}


		}

		
		delete [] state;
		delete [] INDX;
		state=NULL;
		INDX=NULL;



	}
	else
	{
		throw std::logic_error("Set Control proportion over 1");
	}
	
}

void RemoveTrees(int irow, int jcol, double pHealth, double pEinf, double pInf, double pSymp)
{
	double H_trees=Grid::GetCrop(irow, jcol, 0); //Healthy trees
	double E_trees=Grid::GetCrop(irow, jcol, 1); // Infected not infectious
	double I_trees=Grid::GetCrop(irow, jcol, 2); // Infected not symptomatic
	double S_trees=Grid::GetCrop(irow, jcol, 3); // Symptomatic
	
	H_trees=H_trees*(1-pHealth); //healthy
	E_trees=E_trees*(1-pEinf); //infected 
	I_trees=I_trees*(1-pInf); //infectious
	S_trees=S_trees*(1-pSymp); //Symp
	double R_trees = 1 - H_trees - E_trees-I_trees - S_trees;

	Grid::SetCrop(irow, jcol, 0, H_trees);
	Grid::SetCrop(irow, jcol, 1, E_trees);
	Grid::SetCrop(irow, jcol, 2, I_trees);
	Grid::SetCrop(irow, jcol, 3, S_trees);
	Grid::SetCrop(irow, jcol, 4, R_trees);
	
	
}


void Distribute2(int irow, int jcol,  int health, double num)
{//distribute the Bugs according to beta with parameter gamma and nu
	if (num == 0) //if there are no Bugs don't worry
	{
		return;
	}
	else
	{
		double P_Sum = 0; // to keep track of probability used so that we readjust
		int numSum = 0;
		std::vector<double> MyDist;
		//copy the distribution vector to MyDist
		FlightP::GetBugFlight(MyDist);

		int NRows = Grid::GetNumRows();
		int NCols = Grid::GetNumCols();
		int Len = MyDist.size();


		//central point 

		double NewV(0);
		double myCrop = Grid::GetCrop(irow, jcol, Grid::GetMaxNumStatus() - 1); //get proportion 'removed'

		int N = 1;
		int M = num;
		double P = MyDist[0];
		int const LR = 1;
		double R[LR];
		int X[1];
		int ifail = 1;
		int MODE = 3;
		int const genid(1);
		int const subid(0);
		int lstate(0);
		int *state;
		state = new int[1];
		G05KGF(genid, subid, state, lstate, ifail);
		delete[] state;
		state = new int[lstate];
		ifail = 1;
		G05KGF(genid, subid, state, lstate, ifail);
		G05TAF(MODE, N, M, P, R, LR, state, X, ifail);
		P_Sum = P_Sum + MyDist[0];
		M = M - X[0];


		if (myCrop<0.9999) //Trees to land on
		{

			NewV = Grid::GetTempBug(irow, jcol, health) + X[0];
			Grid::SetTempBug(irow, jcol, health, NewV);
			numSum = numSum + X[0];
		}
		else
		{
			flyAgain(irow, jcol);
			NewV = Grid::GetTempBug(irow, jcol, health) + X[0];
			Grid::SetTempBug(irow, jcol, health, NewV);
			numSum = numSum + X[0];
		}

		int sideCount(2);
		int LenCount(3);
		while (Len>LenCount - 1) //next square%
		{
			//sides
			int nrow = irow;
			int ncol = jcol + sideCount - 1;
			Reflect(nrow, ncol);

			double myCrop = Grid::GetCrop(nrow, ncol, Grid::GetMaxNumStatus() - 1); //get proportion 'removed'

			if (P_Sum < 1)
			{
				P = (MyDist[LenCount - sideCount]) / (1 - P_Sum);
			}
			else
			{
				P = 1;
			}
			if (P>1)
				P = 1;
			G05TAF(MODE, N, M, P, R, LR, state, X, ifail);
			P_Sum = P_Sum + MyDist[LenCount - sideCount];
			M = M - X[0];

			if (X[0] > 0)
			{
				if (myCrop<0.9999) //Trees to land on
				{
					NewV = Grid::GetTempBug(nrow, ncol, health) + X[0];
					Grid::SetTempBug(nrow, ncol, health, NewV);
					numSum = numSum + X[0];
				}
				else
				{
					flyAgain(nrow, ncol);
					NewV = Grid::GetTempBug(nrow, ncol, health) + X[0];
					Grid::SetTempBug(nrow, ncol, health, NewV);
					numSum = numSum + X[0];
				}
			}



			nrow = irow;
			ncol = jcol - sideCount + 1;
			Reflect(nrow, ncol);
			myCrop = Grid::GetCrop(nrow, ncol, Grid::GetMaxNumStatus() - 1); //get proportion 'removed'

			if (P_Sum < 1)
			{
				P = (MyDist[LenCount - sideCount]) / (1 - P_Sum);
			}
			else
			{
				P = 1;
			}
			if (P>1)
				P = 1;
			G05TAF(MODE, N, M, P, R, LR, state, X, ifail);
			P_Sum = P_Sum + MyDist[LenCount - sideCount];
			M = M - X[0];

			if (X[0] > 0)
			{
				if (myCrop < 0.9999) //Trees to land on
				{
					NewV = Grid::GetTempBug(nrow, ncol, health) + X[0];
					Grid::SetTempBug(nrow, ncol, health, NewV);
					numSum = numSum + X[0];
				}
				else
				{
					flyAgain(nrow, ncol);
					NewV = Grid::GetTempBug(nrow, ncol, health) + X[0];
					Grid::SetTempBug(nrow, ncol, health, NewV);
					numSum = numSum + X[0];
				}
			}

			nrow = irow + sideCount - 1;
			ncol = jcol;
			Reflect(nrow, ncol);
			myCrop = Grid::GetCrop(nrow, ncol, Grid::GetMaxNumStatus() - 1); //Get number in removed section

			if (P_Sum < 1)
			{
				P = (MyDist[LenCount - sideCount]) / (1 - P_Sum);
			}
			else
			{
				P = 1;
			}
			if (P>1)
				P = 1;
			G05TAF(MODE, N, M, P, R, LR, state, X, ifail);
			P_Sum = P_Sum + MyDist[LenCount - sideCount];
			M = M - X[0];

			if (X[0] > 0)
			{
				if (myCrop < 0.9999) //Then there are trees to land on
				{
					NewV = Grid::GetTempBug(nrow, ncol, health) + X[0];
					Grid::SetTempBug(nrow, ncol, health, NewV);
					numSum = numSum + X[0];
				}
				else
				{
					flyAgain(nrow, ncol);
					NewV = Grid::GetTempBug(nrow, ncol, health) + X[0];
					Grid::SetTempBug(nrow, ncol, health, NewV);
					numSum = numSum + X[0];
				}
			}

			nrow = irow - sideCount + 1;
			ncol = jcol;
			Reflect(nrow, ncol);
			myCrop = Grid::GetCrop(nrow, ncol, Grid::GetMaxNumStatus() - 1); //get proportion 'removed'

			if (P_Sum < 1)
			{
				P = (MyDist[LenCount - sideCount]) / (1 - P_Sum);
			}
			else
			{
				P = 1;
			}
			if (P>1)
				P = 1;
			G05TAF(MODE, N, M, P, R, LR, state, X, ifail);
			P_Sum = P_Sum + MyDist[LenCount - sideCount];
			M = M - X[0];

			if (X[0] > 0)
			{
				if (myCrop < 0.9999) //Trees to land on
				{
					NewV = Grid::GetTempBug(nrow, ncol, health) + X[0];
					Grid::SetTempBug(nrow, ncol, health, NewV);
					numSum = numSum + X[0];
				}
				else
				{
					flyAgain(nrow, ncol);
					NewV = Grid::GetTempBug(nrow, ncol, health) + X[0];
					Grid::SetTempBug(nrow, ncol, health, NewV);
					numSum = numSum + X[0];
				}
			}

			//corners

			nrow = irow + sideCount - 1;
			ncol = jcol + sideCount - 1;
			Reflect(nrow, ncol);
			myCrop = Grid::GetCrop(nrow, ncol, Grid::GetMaxNumStatus() - 1); //get proportion 'removed'

			if (P_Sum < 1)
			{
				P = (MyDist[LenCount - 1]) / (1 - P_Sum);
			}
			else
			{
				P = 1;
			}
			if (P>1)
				P = 1;
			G05TAF(MODE, N, M, P, R, LR, state, X, ifail);
			P_Sum = P_Sum + MyDist[LenCount - 1];
			M = M - X[0];
			if (X[0] > 0)
			{
				if (myCrop < 0.9999) //Trees to land on
				{
					NewV = Grid::GetTempBug(nrow, ncol, health) + X[0]; //alice check
					Grid::SetTempBug(nrow, ncol, health, NewV);
					numSum = numSum + X[0];
				}
				else
				{
					flyAgain(nrow, ncol);
					NewV = Grid::GetTempBug(nrow, ncol, health) + X[0];
					Grid::SetTempBug(nrow, ncol, health, NewV);
					numSum = numSum + X[0];
				}
			}
			//NewV=Grid::GetTempBug(nrow, ncol, health)+MyDist[LenCount-1]*num;
			//Grid::SetTempBug(nrow, ncol, health, NewV);

			nrow = irow + sideCount - 1;
			ncol = jcol - sideCount + 1;
			Reflect(nrow, ncol);
			myCrop = Grid::GetCrop(nrow, ncol, Grid::GetMaxNumStatus() - 1); //get proportion 'removed'
			if (X[0] > 0)
			{
				if (P_Sum < 1)
				{
					P = (MyDist[LenCount - 1]) / (1 - P_Sum);
				}
				else
				{
					P = 1;
				}
			}
			if (P>1)
				P = 1;
			G05TAF(MODE, N, M, P, R, LR, state, X, ifail);
			P_Sum = P_Sum + MyDist[LenCount - 1];
			M = M - X[0];
			if (X[0] > 0)
			{
				if (myCrop < 0.9999) //Trees to land on
				{
					NewV = Grid::GetTempBug(nrow, ncol, health) + X[0];
					Grid::SetTempBug(nrow, ncol, health, NewV);
					numSum = numSum + X[0];
				}
				else
				{
					flyAgain(nrow, ncol);
					NewV = Grid::GetTempBug(nrow, ncol, health) + X[0];
					Grid::SetTempBug(nrow, ncol, health, NewV);
					numSum = numSum + X[0];
				}
			}

			nrow = irow - sideCount + 1;
			ncol = jcol + sideCount - 1;
			Reflect(nrow, ncol);
			myCrop = Grid::GetCrop(nrow, ncol, Grid::GetMaxNumStatus() - 1); //get proportion 'removed'

			if (P_Sum < 1)
			{
				P = (MyDist[LenCount - 1]) / (1 - P_Sum);
			}
			else
			{
				P = 1;
			}
			if (P>1)
				P = 1;
			G05TAF(MODE, N, M, P, R, LR, state, X, ifail);
			P_Sum = P_Sum + MyDist[LenCount - 1];
			M = M - X[0];
			if (X[0] > 0)
			{
				if (myCrop < 0.9999) //Trees to land on
				{
					NewV = Grid::GetTempBug(nrow, ncol, health) + X[0];
					Grid::SetTempBug(nrow, ncol, health, NewV);
					numSum = numSum + X[0];
				}
				else
				{
					flyAgain(nrow, ncol);
					NewV = Grid::GetTempBug(nrow, ncol, health) + X[0];
					Grid::SetTempBug(nrow, ncol, health, NewV);
					numSum = numSum + X[0];
				}
			}

			nrow = irow - sideCount + 1;
			ncol = jcol - sideCount + 1;
			Reflect(nrow, ncol);
			myCrop = Grid::GetCrop(nrow, ncol, Grid::GetMaxNumStatus() - 1); //get proportion 'removed'

			if (P_Sum < 1)
			{
				P = (MyDist[LenCount - 1]) / (1 - P_Sum);
			}
			else
			{
				P = 1;
			}
			if (P>1)
				P = 1;
			G05TAF(MODE, N, M, P, R, LR, state, X, ifail);
			P_Sum = P_Sum + MyDist[LenCount - 1];
			M = M - X[0];
			if (X[0] > 0)
			{
				if (myCrop < 0.9999) //Trees to land on
				{
					NewV = Grid::GetTempBug(nrow, ncol, health) + X[0];
					Grid::SetTempBug(nrow, ncol, health, NewV);
					numSum = numSum + X[0];
				}
				else
				{
					flyAgain(nrow, ncol);
					NewV = Grid::GetTempBug(nrow, ncol, health) + X[0];
					Grid::SetTempBug(nrow, ncol, health, NewV);
					numSum = numSum + X[0];
				}
			}
			//diagonalsides

			int Mid = sideCount - 2;

			for (int Mcount = 0; Mcount<Mid; Mcount++)
			{
				nrow = irow + Mcount + 1;
				ncol = jcol + sideCount - 1;
				Reflect(nrow, ncol);
				myCrop = Grid::GetCrop(nrow, ncol, Grid::GetMaxNumStatus() - 1); //get proportion 'removed'

				if (P_Sum < 1)
				{
					P = (MyDist[LenCount - sideCount + Mcount + 1]) / (1 - P_Sum);
				}
				else
				{
					P = 1;
				}

				if (P>1)
					P = 1;

				G05TAF(MODE, N, M, P, R, LR, state, X, ifail);
				P_Sum = P_Sum + MyDist[LenCount - sideCount + Mcount + 1];
				M = M - X[0];
				if (X[0] > 0)
				{
					if (myCrop < 0.9999) //Trees to land on
					{
						NewV = Grid::GetTempBug(nrow, ncol, health) + X[0];
						Grid::SetTempBug(nrow, ncol, health, NewV);
						numSum = numSum + X[0];
					}
					else
					{
						flyAgain(nrow, ncol);
						NewV = Grid::GetTempBug(nrow, ncol, health) + X[0];
						Grid::SetTempBug(nrow, ncol, health, NewV);
						numSum = numSum + X[0];
					}
				}

				nrow = irow - Mcount - 1;
				ncol = jcol + sideCount - 1;
				Reflect(nrow, ncol);
				double myCrop = Grid::GetCrop(nrow, ncol, Grid::GetMaxNumStatus() - 1); //get proportion 'removed'

				if (P_Sum < 1)
				{
					P = (MyDist[LenCount - sideCount + Mcount + 1]) / (1 - P_Sum);
				}
				else
				{
					P = 1;
				}
				if (P>1)
					P = 1;

				G05TAF(MODE, N, M, P, R, LR, state, X, ifail);
				P_Sum = P_Sum + MyDist[LenCount - sideCount + Mcount + 1];
				M = M - X[0];
				if (X[0] > 0)
				{
					if (myCrop < 0.9999) //Trees to land on
					{
						NewV = Grid::GetTempBug(nrow, ncol, health) + X[0];
						Grid::SetTempBug(nrow, ncol, health, NewV);
						numSum = numSum + X[0];
					}
					else
					{
						flyAgain(nrow, ncol);
						NewV = Grid::GetTempBug(nrow, ncol, health) + X[0];
						Grid::SetTempBug(nrow, ncol, health, NewV);
						numSum = numSum + X[0];
					}
				}

				nrow = irow + Mcount + 1;
				ncol = jcol - sideCount + 1;
				Reflect(nrow, ncol);
				myCrop = Grid::GetCrop(nrow, ncol, Grid::GetMaxNumStatus() - 1); //get proportion 'removed'

				if (P_Sum < 1)
				{
					P = (MyDist[LenCount - sideCount + Mcount + 1]) / (1 - P_Sum);
				}
				else
				{
					P = 1;
				}
				if (P>1)
					P = 1;

				G05TAF(MODE, N, M, P, R, LR, state, X, ifail);
				P_Sum = P_Sum + MyDist[LenCount - sideCount + Mcount + 1];
				M = M - X[0];
				if (X[0] > 0)
				{
					if (myCrop < 0.9999) //Trees to land on
					{
						NewV = Grid::GetTempBug(nrow, ncol, health) + X[0];
						Grid::SetTempBug(nrow, ncol, health, NewV);
						numSum = numSum + X[0];
					}
					else
					{
						flyAgain(nrow, ncol);
						NewV = Grid::GetTempBug(nrow, ncol, health) + X[0];
						Grid::SetTempBug(nrow, ncol, health, NewV);
						numSum = numSum + X[0];
					}
				}

				nrow = irow - Mcount - 1;
				ncol = jcol - sideCount + 1;
				Reflect(nrow, ncol);
				myCrop = Grid::GetCrop(nrow, ncol, Grid::GetMaxNumStatus() - 1); //get proportion 'removed'

				if (P_Sum < 1)
				{
					P = (MyDist[LenCount - sideCount + Mcount + 1]) / (1 - P_Sum);
				}
				else
				{
					P = 1;
				}
				if (P>1)
					P = 1;

				G05TAF(MODE, N, M, P, R, LR, state, X, ifail);
				P_Sum = P_Sum + MyDist[LenCount - sideCount + Mcount + 1];
				M = M - X[0];
				if (X[0] > 0)
				{
					if (myCrop < 0.9999) //Trees to land on
					{
						NewV = Grid::GetTempBug(nrow, ncol, health) + X[0];
						Grid::SetTempBug(nrow, ncol, health, NewV);
						numSum = numSum + X[0];
					}
					else
					{
						flyAgain(nrow, ncol);
						NewV = Grid::GetTempBug(nrow, ncol, health) + X[0];
						Grid::SetTempBug(nrow, ncol, health, NewV);
						numSum = numSum + X[0];
					}
				}

				////
				nrow = irow + sideCount - 1;
				ncol = jcol + Mcount + 1;
				Reflect(nrow, ncol);
				myCrop = Grid::GetCrop(nrow, ncol, Grid::GetMaxNumStatus() - 1); //get proportion 'removed'

				if (P_Sum < 1)
				{
					P = (MyDist[LenCount - sideCount + Mcount + 1]) / (1 - P_Sum);
				}
				else
				{
					P = 1;
				}
				if (P>1)
					P = 1;

				G05TAF(MODE, N, M, P, R, LR, state, X, ifail);
				P_Sum = P_Sum + MyDist[LenCount - sideCount + Mcount + 1];
				M = M - X[0];
				if (X[0] > 0)
				{
					if (myCrop < 0.9999) //Trees to land on
					{
						NewV = Grid::GetTempBug(nrow, ncol, health) + X[0];
						Grid::SetTempBug(nrow, ncol, health, NewV);
						numSum = numSum + X[0];
					}
					else
					{
						flyAgain(nrow, ncol);
						NewV = Grid::GetTempBug(nrow, ncol, health) + X[0];
						Grid::SetTempBug(nrow, ncol, health, NewV);
						numSum = numSum + X[0];
					}
				}

				nrow = irow + sideCount - 1;
				ncol = jcol - Mcount - 1;
				Reflect(nrow, ncol);
				myCrop = Grid::GetCrop(nrow, ncol, Grid::GetMaxNumStatus() - 1); //Get proportion removed

				if (P_Sum < 1)
				{
					P = (MyDist[LenCount - sideCount + Mcount + 1]) / (1 - P_Sum);
				}
				else
				{
					P = 1;
				}
				if (P>1)
					P = 1;

				G05TAF(MODE, N, M, P, R, LR, state, X, ifail);
				P_Sum = P_Sum + MyDist[LenCount - sideCount + Mcount + 1];
				M = M - X[0];
				//NewV=Grid::GetTempBug(nrow, ncol, health)+MyDist[LenCount-sideCount+Mcount+1]*num;
				//Grid::SetTempBug(nrow, ncol, health, NewV);
				if (X[0] > 0)
				{
					if (myCrop < 0.9999) //There is crop to land on
					{
						NewV = Grid::GetTempBug(nrow, ncol, health) + X[0];
						Grid::SetTempBug(nrow, ncol, health, NewV);
						numSum = numSum + X[0];
					}
					else
					{
						flyAgain(nrow, ncol);
						NewV = Grid::GetTempBug(nrow, ncol, health) + X[0];
						Grid::SetTempBug(nrow, ncol, health, NewV);
						numSum = numSum + X[0];
					}
				}

				nrow = irow - sideCount + 1;
				ncol = jcol + Mcount + 1;
				Reflect(nrow, ncol);
				myCrop = Grid::GetCrop(nrow, ncol, Grid::GetMaxNumStatus() - 1); //get proportion 'removed'

				if (P_Sum < 1)
				{
					P = (MyDist[LenCount - sideCount + Mcount + 1]) / (1 - P_Sum);
				}
				else
				{
					P = 1;
				}
				if (P>1)
					P = 1;

				G05TAF(MODE, N, M, P, R, LR, state, X, ifail);
				P_Sum = P_Sum + MyDist[LenCount - sideCount + Mcount + 1];
				M = M - X[0];
				if (X[0] > 0)
				{
					if (myCrop < 0.9999) //Trees to land on
					{
						NewV = Grid::GetTempBug(nrow, ncol, health) + X[0];
						Grid::SetTempBug(nrow, ncol, health, NewV);
						numSum = numSum + X[0];
					}
					else
					{
						flyAgain(nrow, ncol);
						NewV = Grid::GetTempBug(nrow, ncol, health) + X[0];
						Grid::SetTempBug(nrow, ncol, health, NewV);
						numSum = numSum + X[0];
					}
				}

				nrow = irow - sideCount + 1;
				ncol = jcol - Mcount - 1;
				Reflect(nrow, ncol);
				myCrop = Grid::GetCrop(nrow, ncol, Grid::GetMaxNumStatus() - 1); //get proportion 'removed'

				if (P_Sum < 1)
				{
					P = (MyDist[LenCount - sideCount + Mcount + 1]) / (1 - P_Sum);
				}
				else
				{
					P = 1;
				}
				if (P>1)
					P = 1;
				G05TAF(MODE, N, M, P, R, LR, state, X, ifail);
				P_Sum = P_Sum + MyDist[LenCount - sideCount + Mcount + 1];
				M = M - X[0];
				if (X[0] > 0)
				{
					if (myCrop < 0.9999) //Trees to land on
					{
						NewV = Grid::GetTempBug(nrow, ncol, health) + X[0];
						Grid::SetTempBug(nrow, ncol, health, NewV);
						numSum = numSum + X[0];
					}
					else
					{
						flyAgain(nrow, ncol);
						NewV = Grid::GetTempBug(nrow, ncol, health) + X[0];
						Grid::SetTempBug(nrow, ncol, health, NewV);
						numSum = numSum + X[0];
					}
				}


			}


			sideCount++; //these keep track of where you are on square see Fig 3 ECBSpatialModel.pdf
			LenCount = LenCount + sideCount;

		}
		delete[] state;
		state = NULL;
		if (numSum < num)
		{
			NewV = Grid::GetTempBug(irow, jcol, health) + num - numSum;
			Grid::SetTempBug(irow, jcol, health, NewV);
		}
	}

		

		
}

void Reflect(int& nrow, int& ncol)
{
	int NRows=Grid::GetNumRows();
	int NCols=Grid::GetNumCols(); 
	while((nrow>NRows-1)||(nrow<0))
	{
		if (nrow<0)
		{
		nrow=-1-nrow;
		}
		else if (nrow>NRows-1)
		{
		nrow=2*NRows-1-nrow;
		}
	}
	while((ncol>NCols-1)||(ncol<0))
	{
		if (ncol<0)
		{
			ncol=-1-ncol;
		}
		else if
		(ncol>NCols-1)
		{
		ncol=2*NCols-1-ncol;
		}
	}


}

void flyAgain(int& nrow, int& ncol)
{
	//Find a near by field of maize
	
	int maxDis=Grid::GetNumCols();
	if (maxDis<Grid::GetNumRows()) maxDis=Grid::GetNumRows();
	for (int kcount=1; kcount<maxDis; kcount++)
	for (int irow=nrow-kcount; irow<nrow+kcount+1; irow++)
	{
		for (int jcol=ncol-kcount; jcol<ncol+kcount+1; jcol++)
		{
			if ((irow<Grid::GetNumRows())&&(irow>-1)&&(jcol<Grid::GetNumCols())&&(jcol>-1))
			{
				if (Grid::GetCrop(irow, jcol, Grid::GetMaxNumStatus()-1)<0.999999)
				{
					nrow=irow;
					ncol=jcol;
					return;
				}
			}

		}
	}
	
	throw std::logic_error("No Oranges for the Bugs :-("); 


}

void WriteCrops(int iyr)
{
	char path_buffer[_MAX_PATH];
	char Npath_buffer[_MAX_PATH];
	_getcwd( path_buffer, _MAX_PATH );
	strcpy_s(Npath_buffer, path_buffer);
	char Crop[50];
    int n=sprintf_s(Crop, 50, "\\OutFiles\\TreeH%d.txt", iyr);
	strcat_s(Npath_buffer, Crop);
	std::ofstream OutF(Npath_buffer);


	strcpy_s(Npath_buffer, path_buffer);
	n=sprintf_s(Crop, 50, "\\OutFiles\\TreeE%d.txt", iyr); //infected not infectious
	strcat_s(Npath_buffer, Crop);
	std::ofstream OutFE(Npath_buffer);
	
	
	strcpy_s(Npath_buffer, path_buffer);
	n=sprintf_s(Crop, 50, "\\OutFiles\\TreeI%d.txt", iyr);
	strcat_s(Npath_buffer, Crop);
	std::ofstream OutFI(Npath_buffer);

	strcpy_s(Npath_buffer, path_buffer);
	n=sprintf_s(Crop, 50, "\\OutFiles\\TreeS%d.txt", iyr);
	strcat_s(Npath_buffer, Crop);
	std::ofstream OutFS(Npath_buffer);



	if ((OutF)&&(OutFI)&&(OutFS))
	{
		
			for (int jcount=0; jcount<Grid::GetNumRows(); jcount++)
			{
				for (int icount=0; icount<Grid::GetNumCols(); icount++)
				{
					OutF<<Grid::GetCrop(jcount, icount, 0)<<'\t';
					OutFE<<Grid::GetCrop(jcount, icount, 1)<<'\t';
					OutFI<<Grid::GetCrop(jcount, icount, 2)<<'\t';
					OutFS<<Grid::GetCrop(jcount, icount, 3)<<'\t';

				}
				OutF<<'\n';
				OutFE<<'\n';
				OutFI<<'\n';
				OutFS<<'\n';
			}
			
		
	}
	else
		throw std::logic_error("Error in write crops");

}

void WriteBugs(int imth)
{
	char path_buffer[_MAX_PATH];
	char Npath_buffer[_MAX_PATH];
	_getcwd( path_buffer, _MAX_PATH );
	strcpy_s(Npath_buffer, path_buffer);
	char ECBStr[50];
    int n=sprintf_s(ECBStr, 50, "\\OutFiles\\Bugs%d.txt", imth);
	strcat_s(Npath_buffer, ECBStr);
		
	std::ofstream OutF(Npath_buffer);
	if (OutF)
	{
		for (int jcount=0; jcount<Grid::GetNumRows(); jcount++)
		{
			for (int icount=0; icount<Grid::GetNumCols(); icount++)
			{
		
			    double totHPop, totIPop;
				//Grid::GetOWLav(jcount, icount, imth,totHPop, totIPop);
				totIPop = Grid::GetNumG(jcount, icount, 1); //healthy 
				//OutF<<totHPop+totIPop<<'\t';
				OutF << totIPop << '\t';

			}
			OutF<<'\n';
		}
	}
	else
		throw std::logic_error("Error in write bugs");

}


void WriteBugsT(int imth)
{
	char path_buffer[_MAX_PATH];
	char Npath_buffer[_MAX_PATH];
	_getcwd( path_buffer, _MAX_PATH );
	strcpy_s(Npath_buffer, path_buffer);
	char ECBStr[50];
    int n=sprintf_s(ECBStr, 50, "\\OutFiles\\BugsT%d.txt", imth);
	strcat_s(Npath_buffer, ECBStr);
		
	std::ofstream OutF(Npath_buffer);
	if (OutF)
	{
		for (int jcount=0; jcount<Grid::GetNumRows(); jcount++)
		{
			for (int icount=0; icount<Grid::GetNumCols(); icount++)
			{
		
			    double totHPop, totIPop;
				totHPop=Grid::GetNumG(jcount, icount ,0); //healthy 
				totIPop=Grid::GetNumG(jcount, icount ,1); //healthy 
				OutF<<totHPop+totIPop<<'\t';
				OutF.flush();

			}
			OutF<<'\n';
		}
	}
	else
		throw std::logic_error("Error in write bugs");

}


void WriteControl(int imth)
{
	char path_buffer[_MAX_PATH];
	char Npath_buffer[_MAX_PATH];
	_getcwd( path_buffer, _MAX_PATH );
	strcpy_s(Npath_buffer, path_buffer);
	char ECBStr[50];
    int n=sprintf_s(ECBStr, 50, "\\OutFiles\\Control%d.txt", imth);
	strcat_s(Npath_buffer, ECBStr);
		
	std::ofstream OutF(Npath_buffer);
	if (OutF)
	{
		for (int jcount=0; jcount<Grid::GetNumRows(); jcount++)
		{
			for (int icount=0; icount<Grid::GetNumCols(); icount++)
			{
				int myC=Grid::GetControl(jcount, icount,0);
				OutF<<myC<<'\t';
			}
			OutF<<'\n';
		}
	}
	else
		throw std::logic_error("Error in write bugs");

}

void WriteBelief(int imth)
{
	char path_buffer[_MAX_PATH];
	char Npath_buffer1[_MAX_PATH];
	char Npath_buffer2[_MAX_PATH];
	_getcwd(path_buffer, _MAX_PATH);
	strcpy_s(Npath_buffer1, path_buffer);
    strcat_s(Npath_buffer1, "\\OutFiles\\WorryOfInf.txt");
	strcpy_s(Npath_buffer2, path_buffer);
	strcat_s(Npath_buffer2, "\\OutFiles\\BeliefInControl.txt");

	double tempworry, tempbelief;
	std::ofstream OutF(Npath_buffer1, std::ofstream::app);
	std::ofstream OutF1(Npath_buffer2, std::ofstream::app);
		
	if (OutF1&&OutF)
	{
		int numPlanta = Grid::GetPlantaInfoSize(); //returns the number of groves in the landscape
		for (int icount = 0; icount<numPlanta; icount++)
		{
			Grid::GetPlantaScore(icount, tempworry, tempbelief);   //Get worry status and belief in control
			OutF << tempworry << '\t';
			OutF1 << tempbelief << '\t';
		}
		OutF << '\n';
		OutF1 <<  '\n';
		OutF.flush();
		OutF1.flush();

	}
	else
		throw std::logic_error("Error in write belief");

}


void FarmerChoice(int numInd, double distanceRange, double OpinionRange, int frequControl, int HistControl, double ReductionBelief, double impulse)
//numInd are the number of growers I talk to,  distanceRange distance parameter for distance probabilities, OpinionRange is distance parameter for opinions
//Frequency of control is frequency of control specied by CHMA co-ordinator (number of imes a year)
//HistControl - how far back farmers look at control when assessing whether control program is working
//ReductionBelief the propotion of belief in control lost if on average more than one % of trees become sympotomatic and recomended frequency of control was applied during the last HistControl months
//Impulse is the relative weight to which we pay attention to Expert 
{
	//Get last times level of worry 
	std::vector<double> myWorries;
	std::vector<double> myNewWorries;
	//Get last times level of belief in control 
	std::vector<double> myBelief;
	std::vector<double> myNewBelief;
	//Grower and distance weight
	struct GWeights
	{
		int growerID;
		double Weight;
	};

	double tempworry, tempbelief;
	int numPlanta = Grid::GetPlantaInfoSize(); //returns the number of groves in the landscape
	for (int icount = 0; icount<numPlanta; icount++)
	{
		Grid::GetPlantaScore(icount, tempworry, tempbelief);   //Get worry status and belief in control
		myWorries.push_back(tempworry); // worry about infection
		myBelief.push_back(tempbelief); // belief that control works
	}

	//For each player adjust worry
	for (int pcount = 0; pcount < numPlanta; pcount++)
	{
		
		
			//Find how close each neighbor is to me and select numInd at random based on distance///////////////////////////
			std::vector<GWeights> DisProbs;
			std::vector<int> MyNetwork; // the numInd that you talk to
		
			double sumWeight = 0;
			for (int jcount = 0; jcount < numPlanta; jcount++)
			{
				//work out distance assume location given by first cell
				int myloci = Grid::GetCellLocation(pcount, 0);
				int mylocj = Grid::GetCellLocation(jcount, 0);
				int irow_i, icol_i, irow_j, icol_j;
				GridIDToCoord(myloci, irow_i, icol_i);
				GridIDToCoord(mylocj, irow_j, icol_j);
				double dis = sqrt(double(pow((irow_i - irow_j), 2) + pow((icol_i - icol_j), 2)));
				//work out covariance measure
				//exponential covariance model 
				double r = dis / distanceRange;
				double myweight = exp(-r);
				GWeights tempDis;
				tempDis.growerID = jcount;
				tempDis.Weight = myweight;
				DisProbs.push_back(tempDis);
				sumWeight = sumWeight + myweight;
			}

			//Start to select the people I listen to
			MyNetwork.push_back(pcount); //The first person in my network is me
			sumWeight = sumWeight - DisProbs[pcount].Weight; //adjust the sum of distances to omit me (should be no different)
			DisProbs.erase(DisProbs.begin() + pcount);

		
			//randomly select numInd from a uniform distribution
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
			double *urand;
			int NumVals = numInd;
			urand = new double[NumVals];
			ifail = 1;
			G05SAF(NumVals, state, urand, ifail);

			if (sumWeight < 0)
				throw std::logic_error("Error in farmer choice at 1");

			//This section samples numInd players without replacement according to the PDF described by the distance weighting
			for (int iicount = 0; iicount < numInd; iicount++)
			{
				double myPos = sumWeight*urand[iicount];  //defines point in empirical distribution from whcih we get player ID

				//Where is that in our distribution of weights
				double TempsumWeight(0);

				if (DisProbs.size() < 1)
					throw std::logic_error("Error in farmer choice at 2");

				int CheckPerSonAdded(0);

				for (int jcount = 0; jcount < DisProbs.size(); jcount++)
				{
					double UpLim = TempsumWeight + DisProbs[jcount].Weight;
					if ((myPos>TempsumWeight) && (myPos <= TempsumWeight + UpLim)) //Then the random selection picks this player
					{
						int myPlayer = DisProbs[jcount].growerID;
						MyNetwork.push_back(myPlayer);
						sumWeight = sumWeight - DisProbs[jcount].Weight;
						DisProbs.erase(DisProbs.begin() + jcount);
						CheckPerSonAdded = 1;
						jcount = DisProbs.size() + 1; //jump out of loop 
					}
					else
					{
						TempsumWeight = UpLim;
					}
				}
				if (CheckPerSonAdded<0.5)
					throw std::logic_error("Error in farmer choice at 3");
			} //iicount

			delete urand;
			urand = NULL;
			//MyNetwork should now be numInd +1 playes long as it includes me


			//Update worry based on closeness of opinion

			std::vector<double> tempWeights;
			sumWeight = 0;
			double MyConfidence = 0.25; // this value is between 0 and 1 and defines how much a grower maintains thier own opinion. 
			for (int jcount = 0; jcount < numInd + 1; jcount++)
			{
				//work out distance assume location given by first cell
				int ncount = MyNetwork[jcount]; //neighbour number
				double dis = sqrt(pow(myWorries[pcount] - myWorries[ncount], 2));
				//work out covariance measure
				//exponential covariance model 
				//double r = dis / 0.01;
				double r = dis / OpinionRange;
				double myweight = exp(-r);
				tempWeights.push_back(myweight);
				sumWeight = sumWeight + myweight;
			}
			//standardise weights and adjust opinions
			double newOpinion = MyConfidence * myWorries[pcount] * (1 - impulse);
			//Add in the Expert factor here. 
			//Expert has a relative weighting of 'impulse' so of total weights this is sumWeight*impulse

			if (impulse > 0)
			{
				double Imprand[1];
				double var = 0.01;
				G05SKF(1, 0, var, state, Imprand, ifail);
				double tmp_impulse = impulse + Imprand[0];
				if (tmp_impulse < 0)
					tmp_impulse = 0;
				//	sumWeight = sumWeight + impulse*sumWeight;
				for (int jcount = 0; jcount < numInd + 1; jcount++)
				{
					tempWeights[jcount] = (1 - MyConfidence)*(1 - tmp_impulse)*tempWeights[jcount] / sumWeight;
					int ncount = MyNetwork[jcount]; //neighbour number
					newOpinion = newOpinion + tempWeights[jcount] * myWorries[ncount];
				}
				newOpinion = newOpinion + tmp_impulse * 1; //This is weight of expert (i.e. impulse)*worry which is 1

			}
			else
			{
				for (int jcount = 0; jcount < numInd + 1; jcount++)
				{
					tempWeights[jcount] = tempWeights[jcount] / sumWeight * (1 - MyConfidence);
					int ncount = MyNetwork[jcount]; //neighbour number
					newOpinion = newOpinion + tempWeights[jcount] * myWorries[ncount];
				}
			}




			//adjust new opinion again for personal experience. If they have disease on the plantation then worry goes up to 0.9
			int persWorry;
			Grid::GetPlantaDisStatus(pcount, persWorry); //this returns a 1 if the syptoms affect more than 5% of the trees
			if (persWorry > 0)
				newOpinion = 1.0;
			myNewWorries.push_back(newOpinion);
			tempWeights.clear();

			///////////Belief Section///////////////////////

			sumWeight = 0;
			for (int jcount = 0; jcount < numInd + 1; jcount++)
			{
				int ncount = MyNetwork[jcount]; //neighbour number
				double dis = sqrt(pow(myBelief[pcount] - myBelief[ncount], 2));
				//work out covariance measure
				//exponential covariance model 
				double r = dis / OpinionRange;
				double myweight = exp(-r);
				tempWeights.push_back(myweight);
				sumWeight = sumWeight + myweight;
			}
			//standardise weights and adjust opinions
			newOpinion = MyConfidence * myBelief[pcount] * (1 - impulse);

			if (impulse > 0)
			{
				double Imprand[1];
				double var = 0.01;
				G05SKF(1, 0, var, state, Imprand, ifail);
				double tmp_impulse = impulse + Imprand[0];
				if (tmp_impulse < 0)
					tmp_impulse = 0;
				//sumWeight = sumWeight + impulse*sumWeight;
				for (int jcount = 0; jcount < numInd + 1; jcount++)
				{
					tempWeights[jcount] = (1 - MyConfidence)* (1 - tmp_impulse)*tempWeights[jcount] / sumWeight;
					int ncount = MyNetwork[jcount]; //neighbour number
					newOpinion = newOpinion + tempWeights[jcount] * myBelief[ncount];
				}
				newOpinion = newOpinion + tmp_impulse * 1; //This is weight of expert(i.e. impulse)*belief which is 1

			}
			else
			{

				for (int jcount = 0; jcount < numInd + 1; jcount++)
				{
					tempWeights[jcount] = (1 - MyConfidence)*tempWeights[jcount] / sumWeight;
					int ncount = MyNetwork[jcount]; //neighbour number
					newOpinion = newOpinion + tempWeights[jcount] * myBelief[ncount];
				}
			}

			//adjust new opinion again for personal experience. If they have disease on the plantation then worry goes up to 0.9
			int persBelief;
			double HLoss = Grid::GetHealthLoss(pcount); //this returns the total proportions of trees that become symptimatic. It must rescaled by the number of cells
			int numCells = Grid::GetNumCells(pcount);
			HLoss = HLoss / double(numCells); // average increase in sympotomatic trees over grove 
			//Did they control
			int Floc = Grid::GetCellLocation(pcount, 0); // if they controlled they did in all cells so can just look in one
			int irow, icol;
			GridIDToCoord(Floc, irow, icol);
			int myCont = 0;
			for (int mcount = 0; mcount < HistControl; mcount++)
			{
				myCont = myCont + Grid::GetControl(irow, icol, mcount); // counts the number of times we have controlled in the last HistControl months
			}

			if (myCont > 0)
				int junk = 0;

			double myRecentCont = double(myCont) / double(HistControl) * 12.0; // this gives an annual frequency based on last HistControl months

			if ((HLoss > 0.01) && (myRecentCont >= frequControl)) //If they controled and loss healthy plants
				newOpinion = ReductionBelief * newOpinion;
			myNewBelief.push_back(newOpinion);
			tempWeights.clear();

			/////////////////////////////////

			DisProbs.clear();
			MyNetwork.clear(); // the numInd that you talk to
			delete[] state;
			state = NULL;


	} //pcount
	



	//randomly select numInd from a uniform distribution
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
	double *urand;
	int NumVals = 2;
	urand = new double[NumVals];
	ifail = 1;
	

	//Choose control based on worry level and reset worry for next time.
	for (int pcount = 0; pcount<numPlanta; pcount++)
	{
		Grid::SetPlantaScore(pcount, myNewWorries[pcount]);// not pcount is index not plantaID which is pcount+1
		Grid::SetPlantaControlBelief(pcount, myNewBelief[pcount]);// not pcount is index not plantaID which is pcount+1

		
		if ((myNewWorries[pcount]>0.6) && (myNewBelief[pcount]>0.4)) 
		{
			//spray plantaion
			int numCells = Grid::GetNumCells(pcount);
			for (int ccount = 0; ccount<numCells; ccount++)
			{
				int Floc = Grid::GetCellLocation(pcount, ccount);
				int irow, icol;
				GridIDToCoord(Floc, irow, icol);
				Grid::SetControl(irow, icol, 1); //This now indicates whether you are in a control program or not

				//This is new code on tree removal that was not used in the paper

				/*double Sym_trees = Grid::GetCrop(irow, icol, 3); // Symptomatic
				if (Sym_trees>WorryThreshold) // WorryThreshold equates to seeing 
				{
					RemoveTrees(irow, icol, 1.0, 1.0, 1.0, 1.0);
				}*/


			}
		}
		else
		{
			int numCells = Grid::GetNumCells(pcount);
			for (int ccount = 0; ccount<numCells; ccount++)
			{
				int Floc = Grid::GetCellLocation(pcount, ccount);
				int irow, icol;
				GridIDToCoord(Floc, irow, icol);
				Grid::SetControl(irow, icol, 0);
			}
		}
		


		Grid::ResetHealthLoss(pcount); //reset healthloss ready for next year
	}

	myWorries.clear();
	myNewWorries.clear();
	myBelief.clear();
	myNewBelief.clear();


}

void GetPropContandInf(double & PropCnt, double & PropInf)
{

	//Status of trees. This is either Healthy=0, Infected but not infectious = 1, infected infectctious no symtoms=2, Symptoms=3, Removed (which includes the area never planted) =4 
	 
	int TotalTrees(0);
	int TotalHealthyCells(0);
	int TotalControl(0);
	for (int icount = 0; icount < Grid::GetNumRows(); icount++)
	{
		for (int jcount = 0; jcount < Grid::GetNumCols(); jcount++)
		{
			double R=Grid::GetCrop(icount, jcount, 4); //is cell empty
			if (R<0.00000001)
			{
				TotalTrees++;
				R = Grid::GetCrop(icount, jcount, 0); //is cell healthy
				if (R>0.99995)
					TotalHealthyCells++;
				
				int MyC=Grid::GetControl(icount, jcount,0);
				if (MyC > 0.5)
					TotalControl++;

			}
		}
	}

	PropInf = 1.0 - double(TotalHealthyCells) / double(TotalTrees);
	PropCnt = double(TotalControl) / double(TotalTrees);

}
