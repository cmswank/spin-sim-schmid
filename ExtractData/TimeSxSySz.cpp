//g++ -o extractroot TimeSxSySz.cpp `root-config --cflags --glibs`
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h> 
#include <stdlib.h>
#include <math.h>
#include <cmath> 
#include "TFile.h"
#include "TTree.h"
using namespace std;
#define Bin 10001 //this should add one to the par file. 
#define Event 100 

void twoDvec2txt(const char* filename,double sx[][Bin],double sy[][Bin], double sz[][Bin],double t[][Bin],double x[][Bin],double y[][Bin],double z[][Bin]){     
  ofstream  myfile;
  myfile.open (filename);
      myfile <<"Sx"<<" Sy"<<" Sz"<<" x"<<" y"<<" z"<<" t"<<" Event num: "<<Event<<" Bin num: "<<Bin<<"\n";	  
for(int i=0; i<Event; i++){
    //myfile<<"Sx ";
  	for(int j=0;j<Bin;j++){

    myfile << sx[i][j] <<" ";
      //myfile<<"Sy ";
    myfile << sy[i][j] <<" ";
    //myfile<<"Sz ";
    myfile << sz[i][j] <<" ";
	myfile << x[i][j] << " ";
	myfile << y[i][j] << " ";
	myfile << z[i][j] << " ";
	//myfile<<"t\n ";
    myfile<<t[1][j]<<" ";
	myfile<<"\n";
	}
  }
	myfile.close();
  
}

 	
int main()
{	
    const double pi = 3.1415926535897;  //pi
   ///The following code gets the data from root file and exports it to a text file... or maybe I'll export it to a new root file with a branch name thats really hard to spell. 

    TFile *myFile = TFile::Open("/home/cmswank/spin-sim/ExtractData/run9.root");
    TTree *dat = (TTree*)myFile->Get("dat");
    // Create a TTreeReader for the tree, for instance by passing the
    // TTree's name and the TDirectory / TFile it is in.
    //TTreeReader myReader("ntuple", myFile);
    //TTreeReader myReader2("anaTree", myFile2); //TTreeReader doesn't seem to work with nicely with g++ 4.8 
    // Int_t linum=dat->GetEntries();//The number of events?
    Int_t N;
    // printf("linum=%d\n",linum);
    //  int landn=linum*N;
    //  printf("landn=%d\n",landn);
    
    Double_t Sx[Bin];
    Double_t Sy[Bin];     
    Double_t Sz[Bin];
    Double_t t[Bin];
	Double_t X[Bin];
    Double_t Y[Bin];     
    Double_t Z[Bin];




    //  double *Spinx[Bin],*Spiny[Bin],*Spinz[Bin],*tn[Bin];
    double (*Spinx)[Bin]=new double[Event][Bin];
    double (*Spiny)[Bin]=new double[Event][Bin];
    double (*Spinz)[Bin]=new double[Event][Bin];
    double (*tn)[Bin]=new double[Event][Bin];
	double (*posx)[Bin]=new double[Event][Bin];
	double (*posy)[Bin]=new double[Event][Bin];
	double (*posz)[Bin]=new double[Event][Bin];



    for(int ii=0;ii<Event;ii++)
    {
      dat->GetEntry(ii);
      //dat->SetBranchAddress("n",&N);
      // printf("n=%d\n",N);
      for(int j=0;j<Bin;j++){
      	dat->SetBranchAddress("t",&t);
       	dat->SetBranchAddress("sx",&Sx);
       	dat->SetBranchAddress("sy",&Sy);
	    dat->SetBranchAddress("sz",&Sz);
		dat->SetBranchAddress("x",&X);
       	dat->SetBranchAddress("y",&Y);
	   	dat->SetBranchAddress("z",&Z);
	tn[ii][j]=t[j];
	Spinx[ii][j]=Sx[j];
	Spiny[ii][j]=Sy[j];
	Spinz[ii][j]=Sz[j];
	posx[ii][j]=X[j];
	posy[ii][j]=Y[j];
	posz[ii][j]=Z[j];
     }
    } 
    // cout<<ii<<"\n";	
   //Export function.
    twoDvec2txt("/home/cmswank/spin-sim/ExtractData/dataSpecular.txt",Spinx,Spiny,Spinz,tn,posx,posy,posz);
    return 0;
}
