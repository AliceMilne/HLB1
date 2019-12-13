#ifndef __Model_H__
#define __Model_H__

void InitialiseModel(int method, double alpha_W, double beta_W, double alpha_B, double beta_B); //  proportion that believe control works
void SetLand(int readin, double alpha_W, double beta_W, double alpha_B, double beta_B); // 0 to read from file and 1 to create, proportion that believe control works
void DeleteModel();

void ModelMonth(int imth, double randInf, double kill_rate, int frequControl); //if randInf=0 then no random infection event if 1 then a randomly chosed cell has a single infected bug land with probability randIf
void GetPop(int irow, int jcol, double& Bug_H, double& Bug_I);  
void Generation(int irow, int jcol, int iyr, double kill_rate, int frequControl);
void Tree_Generation(int irow, int jcol, int iyr); //progesses tree health status in a time step
void Flight();
void Reflect(int& nrow, int& ncol);
void Distribute2(int irow, int jcol, int health, double num);
void flyAgain(int& nrow, int& ncol);
void SetControl(int imth, double per); //sets control across landscape in month imnth
void RemoveTrees(int irow, int jcol, double pHealth, double infInf, double sInf, double pSymp); // remove a proportion of infected trees from cell irow, jcol


void WriteCrops(int imth); //writes crops for a particular year
void WriteBugs(int imth); //writes infected bugs
void WriteBugsT(int imth); //writes total bugs 
void WriteBelief(int imth);

void WriteControl(int imnth);


void FarmerChoice(int numInd, double distanceRange, double OpinionRange, int frequControl, int HistControl, double ReductionBelief, double impulse);
void GetPropContandInf(double & PropCnt, double & PropInf);



#endif