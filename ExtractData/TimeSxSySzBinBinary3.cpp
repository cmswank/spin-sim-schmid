//g++ -o extractrootbinary TimeSxSySzBinBinary3.cpp `root-config --cflags --glibs`
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h> 
#include <stdlib.h>
#include <math.h>
#include <cmath> 
#include "TFile.h"
#include "TTree.h"
#define mem 1000000
using namespace std;

//#define Bin 1001 //this should add one to the par file. 
//#define Event 1000 

void twoDvec2txt(const char* filename,double sx[][mem],double sy[][mem], double sz[][mem],double t[][mem],double x[][mem],double y[][mem],double z[][mem], double vx[][mem], double vy[][mem], double vz[][mem],const int Bin, const int Event){     
  ofstream  myfile;
   
 myfile.open (filename,ios::binary | ios::out);
     // myfile <<"Sx"<<" Sy"<<" Sz"<<" x"<<" y"<<" z"<<" t"<<" Event num: "<<Event<<" Bin num: "<<Bin<<"\n";	  
for(int i=0; i<Event; i++){
    //myfile<<"Sx ";
  	for(int j=0;j<Bin;j++){
      double Sx; 
      double Sy;
	double Sz;
	double X;
	double Y;
	double Z;
	double VX;
	double VY;
	double VZ;
	double T;

	Sx= sx[i][j];
	
		myfile.write(reinterpret_cast <const char*> (&Sx),sizeof(double));      
//myfile<<"Sy ";
	  Sy= sy[i][j];
	myfile.write(reinterpret_cast <const char*> (&Sy), sizeof(double)); 
    //myfile<<"Sz ";
	   Sz= sz[i][j];
	X= x[i][j];
		Y= y[i][j];
	Z= z[i][j];
	VX=vx[i][j];
	VY=vy[i][j];
	VZ=vz[i][j];

	T=t[1][j];
	myfile.write(reinterpret_cast <const char*> (&Sz),sizeof(double)); 
	myfile.write(reinterpret_cast <const char*> (&X), sizeof(double)); 
	myfile.write(reinterpret_cast <const char*> (&Y), sizeof(double)); 
	myfile.write(reinterpret_cast <const char*> (&Z), sizeof(double)); 
	myfile.write(reinterpret_cast <const char*> (&VX), sizeof(double)); 
	myfile.write(reinterpret_cast <const char*> (&VY), sizeof(double)); 
	myfile.write(reinterpret_cast <const char*> (&VZ), sizeof(double)); 


//	myfile << x[i][j] << " ";
//	myfile << y[i][j] << " ";
//	myfile << z[i][j] << " ";
	myfile.write(reinterpret_cast <const char*> (&T), sizeof(double)); 	
	//myfile<<"t\n ";
    	//myfile<<t[1][j]<<" ";
	//myfile<<"\n";
	}
 }
	myfile.close();
  
}

 	
int main(int argc, char *argv[] )
{	string s[argc];
	int params[argc];
	int params1;
	string s1;	
	if (argc!=3) cout<<"Remember its Run #, Bin #, Event#";
	assert(argc!=3);
	
	for (int i = 0; i<3; i++){
	s[i]=argv[i];
	s1=s[i];
    stringstream convert(s1);//object from the class stringstream.
    convert>>params1; //the object has the value convert and streams it
	params[i]=params1;
	}
	const int run=params[0]; 
    const int Bin=params[1]; 
	const int Event=params[2];
	
    const double pi = 3.1415926535897;  //pi
   ///The following code gets the data from root file and exports it to a text file... or maybe I'll export it to a new root file with a branch name thats really hard to spell. 

    TFile *myFile = TFile::Open("/home/cmswank/spin-sim-heatflush/ExtractData/run"<<argv[0]<<".root");
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
    Double_t VX[Bin];
    Double_t VY[Bin];
    Double_t VZ[Bin];


    //  double *Spinx[Bin],*Spiny[Bin],*Spinz[Bin],*tn[Bin];
    double (*Spinx)[Bin]=new double[Event][Bin];
    double (*Spiny)[Bin]=new double[Event][Bin];
    double (*Spinz)[Bin]=new double[Event][Bin];
    double (*tn)[Bin]=new double[Event][Bin];
	double (*posx)[Bin]=new double[Event][Bin];
	double (*posy)[Bin]=new double[Event][Bin];
	double (*posz)[Bin]=new double[Event][Bin];

	double (*velx)[Bin]=new double[Event][Bin];
	double (*vely)[Bin]=new double[Event][Bin];
	double (*velz)[Bin]=new double[Event][Bin];


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
		dat->SetBranchAddress("vx",&VX);
		dat->SetBranchAddress("vy",&VY);
		dat->SetBranchAddress("vz",&VZ);



		tn[ii][j]=t[j];
		Spinx[ii][j]=Sx[j];
		Spiny[ii][j]=Sy[j];
		Spinz[ii][j]=Sz[j];
		posx[ii][j]=X[j];
		posy[ii][j]=Y[j];
		posz[ii][j]=Z[j];
		velx[ii][j]=VX[j];
		vely[ii][j]=VY[j];
		velz[ii][j]=VZ[j];


     }
    } 
    // cout<<ii<<"\n";	
   //Export function.
   twoDvec2txt("/home/cmswank/spin-sim-heatflush/ExtractData/OutputHeatflush"<<argv[0]<<".dat",Spinx,Spiny,Spinz,tn,posx,posy,posz,velx,vely,velz,Bin,Event);
    return 0;
}
