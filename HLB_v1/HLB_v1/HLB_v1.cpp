// HLB_v1.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <stdexcept>
#include  <fstream>
#include  <iostream>
#include <vector>
#include <direct.h>
#include "HLBModel.h"

int _tmain(int argc, _TCHAR* argv[])
{
	//Initialise the landscape
	try
	{
		int iniMode = 1; //0 = uniform trees no landscape structure; 1= read from files; 2= uniform but has plantations and blocks
		int MaxMonth = 36; // number of months that we model
		int FarmDes = 1; //flag to switch between model run and farmer decision run
		double RandIf = 0.0; // The probability that a single bug arrives at a randomly chosen host cell
		//int MyInfluence(1); //Defines what isets the influence of farmers over each other. if MyInfluence=0 then distance if 1 the closeness of opinion
		//Farmer choice parameters and others that we might change

		int numInd = 60; // number of individuals each grower talks to
		double distanceRange = 30; //range parameter for PDF that describes change of each grower being in the number of individuals spoken to
		double OpinionRange = 0.05; // range parameter for opinion dynamics
		int frequControl = 12; //frequency of control recommneded by control program (currently only works for 12)
		int HistControl = 2; //how far back a grower looks when considering whether control is working
		double ReductionBelief = 0.2; //proportion by which belief is reduced
		double alpha_W = 1;
		double beta_W = 1;
		double alpha_B = 5;
		double beta_B = 7;
		double ExpertImpulse = 0.2; // the weighting that a Time opinion might have
		int ExpertImpulseFrequ = 0; // the frequency with which Expert talks to growers
		double kill_rate = 0.94;
		int offset = 1; 

		char path_buffer[_MAX_PATH];
		char Npath_buffer[_MAX_PATH];
		_getcwd(path_buffer, _MAX_PATH);
		strcpy_s(Npath_buffer, path_buffer);
		char Crop[50];
		strcat_s(Npath_buffer, "\\BigSimulationResults.txt");

		std::ofstream WriteF(Npath_buffer);

		int numInd_v[3] = { 10, 60, 100 };
		double distanceRange_v[3] = { 10, 30, 120 };
		double OpinionRange_v[3] =  { 0.05, 0.2, 0.5 }; //{ 0.01, 0.2, 0.5 };
		int ExpertImpulseFrequ_v[4] = { 0, 3, 6, 12 }; // the frequency with which Expert talks to growers
		double ExpertImpulse_v[3] = { 0.2, 0.4, 0.6 }; // the weighting that a Time opinion might have
		int HistControl_v[3] = { 2, 6 };
		double ReductionBelief_v[3] = { 0.2, 0.4, 0.6 }; //proportion by which belief is reduced
		double mu_W[3] = { 0.2, 0.5, 0.8 };
		double Var_W[3] = { 0.01, 0.05, 0.1 }; //{ 0.01, 0.03, 0.05 };
		double mu_B[3] = { 0.2, 0.5, 0.8 };
		double Var_B[3] = { 0.01, 0.05, 0.1 };  //{ 0.01, 0.03, 0.05 };
		int frequControl_v[3] = { 4, 6, 12 }; // 
		double kill_rate_v[2] = { 0.7, 0.94 }; //{ 0.6, 0.94 };

													
		int icount1 = 1;
		numInd = numInd_v[icount1];
		int icount2 = 1;
		distanceRange = distanceRange_v[icount2];
		int icount3 = 0;  // 0
		OpinionRange = OpinionRange_v[icount3];
        int icount4 = 0; 
		ExpertImpulseFrequ = ExpertImpulseFrequ_v[icount4];
		int icount5 = 0; 
		ExpertImpulse = ExpertImpulse_v[icount5];
		int icount6 = 0; 
		HistControl = HistControl_v[icount6];
		int icount7 = 1;
		ReductionBelief = ReductionBelief_v[icount7]; //This doesn't seem to have a big effect so removed for speed
		int icount8 = 0;
		int icount9 = 2;
		double mu = mu_W[icount8];
		double V = Var_W[icount9];
		alpha_W = (mu*mu*(1 - mu) - V*mu) / V;
		beta_W = (mu*(1 - mu)*(1 - mu) - V*(1 - mu)) / V;
	//	if ((alpha_W < 1) && (beta_W < 1))
		//	throw std::logic_error("Error on Beta distribution parameters");
        int icount10 = 1; 
		int icount11 = 2; 
		mu = mu_B[icount10];
		V = Var_B[icount11];
		alpha_B = (mu*mu*(1 - mu) - V*mu) / V;
		beta_B = (mu*(1 - mu)*(1 - mu) - V*(1 - mu)) / V;
	//	if ((alpha_B < 1) && (beta_B < 1))
		//	throw std::logic_error("Error on Beta distribution parameters");
        int icount12 = 2; 
		frequControl = frequControl_v[icount12];
		int icount13 = 1; 
		kill_rate = kill_rate_v[icount13];

		InitialiseModel(iniMode, alpha_W, beta_W, alpha_B, beta_B);
		WriteCrops(0);
//		char path_buffer[_MAX_PATH];
//		char Npath_buffer[_MAX_PATH];
		_getcwd(path_buffer, _MAX_PATH);
		strcpy_s(Npath_buffer, path_buffer);
		strcat_s(Npath_buffer, "\\OutFiles\\MyErrors.txt");
		std::ofstream of(Npath_buffer);

		char Npath_buffer1[_MAX_PATH];
		char Npath_buffer2[_MAX_PATH];
		strcpy_s(Npath_buffer1, path_buffer);
		strcat_s(Npath_buffer1, "\\OutFiles\\WorryOfInf.txt");
		strcpy_s(Npath_buffer2, path_buffer);
		strcat_s(Npath_buffer2, "\\OutFiles\\BeliefInControl.txt");
		std::ofstream OutF(Npath_buffer1);
		std::ofstream OutF1(Npath_buffer2);

		of.close();

		if (FarmDes == 0)
		{
			//This code runs the model with set proportions of control
			for (int icount = 0; icount < MaxMonth; icount++)
			{
				SetControl(icount, 0.0); //Sets control and might also removes trees and bugs with them one day. if per < 0 reads from file otherwise random allocation of that proportion of spray
				//std::cout<<"ModelStar"<<'\n';
				ModelMonth(icount, RandIf, kill_rate, frequControl);
				std::cout << "ModelDone" << icount << '\n';
				WriteCrops(icount + 1);
				WriteBugs(icount); // this function writes a map of total bugs
				WriteControl(icount);
			}
		}
		else
		{
			//This code runs the model with farmer decisions

			for (int icount = 0; icount<MaxMonth; icount++)
			{
				WriteBelief(icount);
				//			FarmerChoice(MyInfluence); //Sets control and might also removes trees and bugs with them one day. if per < 0 reads from file otherwise random allocation of that proportion of spray

				double impulse = 0;
				if (ExpertImpulseFrequ>0)
				{
					double denom = double(ExpertImpulseFrequ) / 12.0;
					double ans1 = double(icount + offset)*denom;
					ans1 = ans1 - floor(ans1);
					if (ans1 < 0.0000001)
						impulse = ExpertImpulse;
				}

				FarmerChoice(numInd, distanceRange, OpinionRange, frequControl, HistControl, ReductionBelief, impulse);
				ModelMonth(icount, RandIf, kill_rate, frequControl);
				WriteCrops(icount);
				WriteBugs(icount); // this function writes a map of infected bugs for total see Tbugs
				WriteControl(icount);
				std::cout << "DecisionDone" << icount << '\n';
				double PropCnt, PropInf;
				GetPropContandInf(PropCnt, PropInf);
			//	if (PropInf > 0.98)
				//	icount = MaxMonth;

			}
			WriteF << numInd << '\t' << distanceRange << '\t' << OpinionRange <<'\t'<< ExpertImpulseFrequ << '\t' << ExpertImpulse << '\t';
			WriteF << HistControl << '\t' << ReductionBelief << '\t' << mu_W[icount8] << '\t' << Var_W[icount9] << '\t';
			WriteF << mu_B[icount10] << '\t' << Var_B[icount11] <<'\t' << frequControl << '\t' << kill_rate << '\t';
			double PropCnt, PropInf;
			GetPropContandInf(PropCnt, PropInf);
			WriteF << PropCnt << '\t' << PropInf << '\n';
			WriteF.flush();

			


		}
		
														
		WriteF.close();

		std::cout << "Job complete";
	

	}
	catch (std::logic_error &E)
	{
		char Mess[500]; 
		strcpy_s(Mess, E.what());
		std::cout<<Mess;
		
		
	}

		return 0;



	
}

