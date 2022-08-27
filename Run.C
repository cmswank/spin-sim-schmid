#include "Run.h"

#include <time.h>
#include <assert.h>
#include <stdlib.h>
#include <fstream>
#include <signal.h>
//#include <fftw3.h>

#include "TMath.h"

#include "BField/BField.h"
//#include "BField/BFieldInterp.h"  REMBER TO CHANGE THE MAKE FILE IF YOU ADDED A CLASS IN A .C FILE (field .o section)
#include "Vector.h"


/**
 * Constructor
 */

RunData::RunData(Int_t n) : id(0), n(0) {
  Initialize(n);
}
/**
 * Destructor
 */
RunData::~RunData() {
  Reset();
}

/**
 * Destructor helper
 */
void RunData::Reset() {
  if (n) {
    n = 0;

    delete[] bounces;
    delete[] scat;
    delete[] steps;
    delete[] backsteps;
    delete[] time;

    delete[] posx;
    delete[] posy;
    delete[] posz;
    delete[] velx;
    delete[] vely;
    delete[] velz;
    delete[] speed;
    delete[] spinx;
    delete[] spiny;
    delete[] spinz;

    delete[] fieldx;
    delete[] fieldy;
    delete[] fieldz;
  }
}

/**
 * Initialize data storage
 */
void RunData::Initialize(Int_t n) {
  // delete arrays if created
  Reset();

  this->n = n;
  this->bounces = new Int_t[n];
  this->scat = new  Int_t[n];
  this->steps = new  Int_t[n];
  this->backsteps = new Int_t[n];
  this->time = new  Double_t[n];
  this->posx = new  Double_t[n];
  this->posy = new  Double_t[n];
  this->posz = new  Double_t[n];
  this->velx = new Double_t[n];
  this->vely = new Double_t[n];
  this->velz = new Double_t[n];
  this->speed = new Double_t[n];
  this->spinx = new Double_t[n];
  this->spiny = new Double_t[n];
  this->spinz = new Double_t[n];

  this->fieldx = new Double_t[n];
  this->fieldy = new Double_t[n];
  this->fieldz = new Double_t[n];
}





/**
 * Create data branches in TTree
 */
void RunData::BranchTTree(TTree *tree) {
  tree->Branch("id", &id, "id/I");
  tree->Branch("n", &n,"n/I");
  tree->Branch("bounces", bounces,"bounces[n]/I");
  tree->Branch("scat", scat,"scat[n]/I");
  tree->Branch("t", time, "t[n]/D");
  tree->Branch("x", posx, "x[n]/D");
  tree->Branch("y", posy, "y[n]/D");
  tree->Branch("z", posz, "z[n]/D");
  tree->Branch("vx", velx, "vx[n]/D");
  tree->Branch("vy", vely, "vy[n]/D");
  tree->Branch("vz", velz, "vz[n]/D");
  tree->Branch("speed", speed, "speed[n]/D");
  tree->Branch("sx", spinx, "sx[n]/D");
  tree->Branch("sy", spiny, "sy[n]/D");
  tree->Branch("sz", spinz, "sz[n]/D");

  tree->Branch("fieldx", fieldx, "fieldx[n]/D");
  tree->Branch("fieldy", fieldy, "fieldy[n]/D");
  tree->Branch("fieldz", fieldz, "fieldz[n]/D");

  tree->Branch("steps", steps, "steps[n]/I");
  tree->Branch("backsteps", backsteps, "backsteps[n]/I");

  char tmp[100];
  sprintf(tmp, "bx[%d]/D", RunData::store_n_bounces-1);
  tree->Branch("bx", bx, tmp);
  sprintf(tmp, "by[%d]/D", RunData::store_n_bounces-1);
  tree->Branch("by", by, tmp);
  sprintf(tmp, "bz[%d]/D", RunData::store_n_bounces-1);
  tree->Branch("bz", bz, tmp);
  sprintf(tmp, "bkind[%d]/D", RunData::store_n_bounces-1);
  tree->Branch("bkind", bkind, tmp);
}

/**
 * Set data branches addresses in TTree
 */
void RunData::SetTTreeBranches(TTree *tree) {
  tree->SetBranchAddress("id", &id);
  tree->SetBranchAddress("n", &n);
  tree->SetBranchAddress("bounces", bounces);
  tree->SetBranchAddress("scat", scat);
  tree->SetBranchAddress("t", time);
  tree->SetBranchAddress("x", posx);
  tree->SetBranchAddress("y", posy);
  tree->SetBranchAddress("z", posz);
  tree->SetBranchAddress("vx", velx);
  tree->SetBranchAddress("vy", vely);
  tree->SetBranchAddress("vz", velz);
  tree->SetBranchAddress("speed", speed);
  tree->SetBranchAddress("sx", spinx);
  tree->SetBranchAddress("sy", spiny);
  tree->SetBranchAddress("sz", spinz);

  tree->SetBranchAddress("fieldx", fieldx);
  tree->SetBranchAddress("fieldy", fieldy);
  tree->SetBranchAddress("fieldz", fieldz);

  tree->SetBranchAddress("steps", steps);
  tree->SetBranchAddress("backsteps", backsteps);

  tree->SetBranchAddress("bx", bx);
  tree->SetBranchAddress("by", by);
  tree->SetBranchAddress("bz", bz);
  tree->SetBranchAddress("bkind", bkind);
}

/**
 * Constructor
 */
ParticleFactory::ParticleFactory() : field(0), boundary(0), scattering(0) {
  for (int i=0; i < 4; i++) {
    vDistribution[i] = 0;
  }
}

/**
 * Create a new particle according to the
 * factory's parameters and return it
 */
Neutron *ParticleFactory::newParticle() {

  newPosition();
  newVelocity();
  newSpin();
   
  // Create the particle
  Neutron *p;
  if (parameters.simType == 1) {
    
    Helium3 *h = new Helium3(temp_pos, temp_vel, temp_spin);
    if (vDistribution[0]) {
      h->vDistribution = vDistribution[0];
    }
    p = (Neutron*)h;
    p->SetGamma(HE3_GAMMA);
    p->SetEDM(0.0);
  }
  else {
    
    p = new Neutron(temp_pos, temp_vel, temp_spin);
    p->SetGamma(NEUTRON_GAMMA);
  }
  if (parameters.simGamma != 0) {
    p->SetGamma(parameters.simGamma);
  }
  if (parameters.simEDM != 0) {
    p->SetEDM(parameters.simEDM);
    //std::cout<<"setting EDM to "<<parameters.simEDM<<"\n";
    //std::cout<<"actual nEDM is "<<p->GetGammaE()<<"\n";
  }
  p->SetBoundary(boundary);
  p->SetField(field);
  p->SetScatter(scattering);
  if (scattering)
    scattering->next();

  // Numerical method:
  p->kutta = parameters.numerical_method;
  if (parameters.numerical_method == 5)
    p->maxError = parameters.max_error;

  p->maxAngle = parameters.max_angle;
  p->no_normalize = parameters.no_normalize;

  return p;
}

/**
 * Create new position vector
 */
Double_t *ParticleFactory::newPosition() {
  boundary->CreateInside(temp_pos);
  return temp_pos;
}

/**
 * Create new velocity vector
 */
Double_t *ParticleFactory::newVelocity() {
  Double_t speed = parameters.speed;
  //std::cout<<speed;
  if (vDistribution[0] && parameters.simType==1)
    speed = vDistribution[0]->GetRandom();
  else if (parameters.vdistrexpo >= 0.) {
    // pick from values distributed as f(v) ~ v^expo
    speed *= pow(gRandom->Rndm(), 1./(1.+parameters.vdistrexpo));
  }
  if (parameters.geometry == 5) {
    // 1-D geometry
    temp_vel[0] = speed;
    temp_vel[1] = 0.;
    temp_vel[2] = 0.;
  }
  else {
    gRandom->Sphere(temp_vel[0],temp_vel[1],temp_vel[2], speed);
    // Apply 2D filter
    if (parameters.geometry == 0 || parameters.geometry == 3) {
      temp_vel[2] = 0;
      Vector::Normalize(temp_vel, speed);
    }
  }
  
  return temp_vel;
}

/**
 * Create new spin vector
 */
Double_t *ParticleFactory::newSpin() {
  // Set the spin (randomize once in any case, for consistency)
  // this way the trajectories are not affected
  gRandom->Sphere(temp_spin[0],temp_spin[1],temp_spin[2],1.);
  if (parameters.spinSettings[1] > -1.1) {
    // single spin polarization
    temp_spin[0] = parameters.spinSettings[0];
    temp_spin[1] = parameters.spinSettings[1];
    temp_spin[2] = parameters.spinSettings[2];
  }
  else {
    // integer values of spinSettings[0]
    Double_t mag;
    int spin_align = (int)(parameters.spinSettings[0]);
    switch (spin_align) {
    case 1:
      // Spin settings (T1) (spin aligned with B0)
      Vector::Copy(parameters.B0, temp_spin);
      break;
    case 2:
    case 4:
      Double_t temp4[3];
      if (spin_align == 2) {
	// T2 random setup (spin perp to B0)
	Vector::Copy(parameters.B0, temp4);
      }
      else {
	// T2 random setup (spin perp to local B)
	BFieldVars var4(0, temp_pos, temp_vel);
	field->getField(temp4, var4);
      }
      Vector::Normalize(temp4);
      mag = Vector::DotProduct(temp_spin, temp4);
      Vector::Normalize(temp4, -mag);
      Vector::Add(temp_spin, temp4, temp_spin);
      break;
    case 3:
      // T1 start aligned with the actual B field
      BFieldVars var5(0, temp_pos, temp_vel);
      field->getField(temp_spin, var5);
      break;
    }
    Vector::Normalize(temp_spin);
  }
  return temp_spin;
}

/**
 * Constructor
 */
Run::Run(Int_t runID) : version_check(1), factory(0), runID(runID), nCurrent(-1),
			field(0), boundary(0), scattering(0),
			mfp(-1), tau_C(-1), dataTree(0), file(0), randomPars(0) {
  sprintf(this->filename, "run%d.root", runID);
  sprintf(this->name, "run%d", runID);

  for (int i=0; i<4; i++)
    vDistribution[i] = 0;
}

/**
 * Destructor
 */
Run::~Run() {
  if (field) {
    delete field;
  }
  for (int i=0; i<4; i++) {
    if (vDistribution[i]) {
      delete vDistribution[i];
    }
  }
  if (factory) {
    if (factory->parameters.uniform_gradient == parameters.uniform_gradient) {
      factory->parameters.uniform_gradient = 0;
    }
    if (factory->parameters.simple_gradient == parameters.simple_gradient) {
      factory->parameters.simple_gradient = 0;
    }
    if (factory->parameters.simple_quad_gradient == parameters.simple_quad_gradient) {
      factory->parameters.simple_quad_gradient = 0;
    }
    delete factory;
  }
}

/**
 * Reload run that has already started
 */
Int_t Run::loadRun() {

  this->file = new TFile(this->filename, "UPDATE");

  this->file->GetObject("dat", dataTree);
  if (!dataTree) {
    this->file->GetObject(this->name, dataTree);
  }
  if (!dataTree) {
    cout << "ERROR reading data tree." << endl;
    return -1;
  }
  this->file->GetObject("par", parameters.tree);
  parameters.LoadParameters();
  
  if (CheckVersion() != 0)
    return -1;
  
  this->file->GetObject("rnd", randomPars);

  // field formulas
  this->file->GetObject("field_formula_x", parameters.field_formula[0]);
  this->file->GetObject("field_formula_y", parameters.field_formula[1]);
  this->file->GetObject("field_formula_z", parameters.field_formula[2]);


  this->nCurrent = dataTree->GetEntries();

  neutronData.SetTTreeBranches(dataTree);

  return 0;

}

/**
 * Check version of data file
 */
int Run::CheckVersion() {

  int version_cmp = 0;

  if (parameters.version[0] > RUN_VERSION[0]) version_cmp =  3;
  else if (parameters.version[0] < RUN_VERSION[0]) version_cmp =  -3;

  if (parameters.version[1] > RUN_VERSION[1]) version_cmp =  2;
  else if (parameters.version[1] < RUN_VERSION[1]) version_cmp =  -2;

  if (parameters.version[2] > RUN_VERSION[2]) version_cmp =  1;
  else if (parameters.version[2] < RUN_VERSION[2]) version_cmp =  -1;
  
  if (version_cmp != 0) {
    if (version_check)
      cout << "ERROR: Version Mismatch:" << endl;
    else
      cout << "Warning: Version Mismatch:" << endl;
    cout << "\tFileVer: " 
	 << parameters.version[0] << "."
	 << parameters.version[1] << "."
	 << parameters.version[2] << endl;
    cout << "\tProgVer: "
	 << RUN_VERSION[0] << "."
	 << RUN_VERSION[1] << "."
	 << RUN_VERSION[2] << endl;
    if (version_check)
      return version_cmp;
  }


  return 0;

}

/**
 * Initialize run
 */
void Run::InitializeRun() {

  neutronData.Initialize(parameters.nBins + 1);

  if (!dataTree) {
    dataTree = createTree(name);
  }

  Double_t gamma = 0;
  switch (parameters.simType) {
  case 0:
    cout << "Neutron Simulation" << endl;
    gamma = NEUTRON_GAMMA;
    break;
  case 1:
    cout << "He3 Simulation" << endl;
    gamma = HE3_GAMMA;
    useBoltzmannDistribution(parameters.temperature);
    break;
  }
  if (parameters.simGamma != 0) {
    gamma = parameters.simGamma;
    cout << "User defined gamma = " << gamma << endl;
  }
  cout.precision(10);
  cout << "gamma = " << gamma << endl;

  // make sure that offset time is positive and < totalTime
  if (parameters.offsetTime < 0. || parameters.offsetTime >= parameters.totalTime) {
    cout << "# invalid offset time provided" << endl;
    parameters.offsetTime = 0.;
  }

  // override the standard gRandom global variable
  // using TRandom3 which is a better generator
  // the generator is seeded with 0 so that it chooses its's 
  // seed based on time??
  delete gRandom;
  if (this->randomPars) {
    gRandom = randomPars;
    cout << "rnd restored" << endl;
  }
  else {
    if (parameters.seed == 0) {
      time_t curtime;
      time(&curtime);
      parameters.seed = (Int_t)curtime;
    }
    gRandom = new TRandom3(parameters.seed);
    parameters.seed = gRandom->GetSeed();
    cout << "seed : " << parameters.seed << endl;
    randomPars = gRandom;
  }

  // Set Geometry
  switch (parameters.geometry) {
  case 0:
    boundary = new CircleBoundary(parameters.radius); // 2D circle
    break;
  case 1:
    // box
    boundary = new BoxBoundary(parameters.boxLow, parameters.boxHigh);
    break;
  case 2:
    boundary = new SphereBoundary(parameters.radius); // sphere
    break;
  case 3:
    boundary = new Box2DBoundary(parameters.boxLow, parameters.boxHigh);
    break;
  case 4:
    boundary = new CylindricalBoundary(parameters.radius, parameters.boxLow[2], parameters.boxHigh[2]);
    break;
  case 5:
    boundary = new Boundary1D(parameters.boxLow[0], parameters.boxHigh[0]);
    break;
  default:
    boundary = 0;
  }

  // Set Uniform Fields
  this->field = new ParticleField(parameters.B0, parameters.E0);
  if (parameters.uniform_gradient)
    this->field->SetUniformGrad(parameters.uniform_gradient);
  else {
    if (parameters.simple_gradient) {
      field->SetSimpleGrad(parameters.simple_gradient);
    }
    if (parameters.simple_quad_gradient) {
      field->SetSimpleQuadGrad(parameters.simple_quad_gradient);
    }
  }

  // Set additional fields
// Set Interpolation fields (generated by COMSOL or other FEM)
///added by C.Swank, TODO include time dependence
if (parameters.Interpparam[0]){
  
  //time independent? 
  if ((int)parameters.Interpparam[0]==1){ 
	  
	double star[3]={(double)parameters.Interpparam[1],(double)parameters.Interpparam[2],(double)parameters.Interpparam[3]};
	double ste[3]={(double)parameters.Interpparam[4],(double)parameters.Interpparam[5],(double)parameters.Interpparam[6]};
	int Ninterp[3]={(int)parameters.Interpparam[7],(int)parameters.Interpparam[8],(int)parameters.Interpparam[9]};
	double nuttin[1]={0.};
    
    std::cout<<"B Field Interpolation"<<endl;
   //constructor for reference. 
   //BFieldInterp::BFieldInterp(bool fromfile, int dir,double *star, double *ste, int N[3], double *datainput)
   //std::cout<<parameters.InterpDir<<"\n";
  	
	if (strcmp(parameters.InterpDir.c_str(),"x")==0)
		{
			this->field->BFIx=new BFieldInterp(true,0,star,ste,Ninterp,nuttin);  
		    this->field->add_field_interp('x', this->field->BFIx);
			std::cout<<"adding Bx interpolation field\n";
			//double tempos[3]={0.3280196619,0.1629835684,-0.8445006203};
			//std::cout<<"test fied "<<this->field->BFIx->interp3D(tempos)<<"\n";
		}
	if (strcmp(parameters.InterpDir.c_str(),"y")==0)
		{
		this->field->BFIy=new BFieldInterp(true,1,star,ste,Ninterp,nuttin);  
		this->field->add_field_interp('y', this->field->BFIy);
		std::cout<<"adding By interpolation field\n";
		}
	if (strcmp(parameters.InterpDir.c_str(),"z")==0)
	{
		this->field->BFIz=new BFieldInterp(true,2,star,ste,Ninterp,nuttin); 
		this->field->add_field_interp('z', this->field->BFIz);
		std::cout<<"adding Bz interpolation field\n";
	}
  } 
}





	if (parameters.SDparam[0]){
	
	if ((int)parameters.SDparam[0]==1){	
		this->field->BSDF=new BDressingFactor();	
		std::cout<<"Spin Dressing active without Modulation"<<endl;
		
	}	
	else if((int)parameters.SDparam[0]==2){
		BDressingCosModFactor *BDCMF= new BDressingCosModFactor();	
		this->field->BSDF=dynamic_cast<BDressingCosModFactor*>(BDCMF);
		this->field->BSDF->fm=parameters.SDparam[3];
    this->field->BSDF->wrf=parameters.SDparam[1];
  		this->field->BSDF->wrf_amp=parameters.SDparam[4];
      this->field->BSDF->t1=parameters.SDparam[5];
      this->field->BSDF->t2=parameters.SDparam[6];
      this->field->BSDF->nmod=parameters.SDparam[7];
      this->field->BSDF->findFM();
      
		std::cout<<"Spin Dressing active with Constant Power Modulation"<<endl;
		
	}
	else if((int)parameters.SDparam[0]==3){
		BDressingCosBModFactor *BDCBMF= new BDressingCosBModFactor();	
		this->field->BSDF=dynamic_cast<BDressingCosBModFactor*>(BDCBMF);
		this->field->BSDF->fm=parameters.SDparam[3];
		this->field->BSDF->amp=parameters.SDparam[4];
	    std::cout<<"Spin Dressing active with Cosine Field Magnitude Modulation"<<endl;
  
	}
	else if((int)parameters.SDparam[0]==4){
		BDressingPulsedBModFactor* BDPBMF = new BDressingPulsedBModFactor();
		this->field->BSDF=dynamic_cast<BDressingPulsedBModFactor*>(BDPBMF);
		this->field->BSDF->fm=parameters.SDparam[3];
  		this->field->BSDF->scale1=parameters.SDparam[4];
  		this->field->BSDF->scale2=parameters.SDparam[5];
  		this->field->BSDF->deltat=parameters.SDparam[6];
		std::cout<<"Spin Dressing active with Pulsed Field Magnitude Modulation"<<endl;
	}
	else if((int)parameters.SDparam[0]==5){
		BDressingPulsedFreqModFactor* BDPFMF = new BDressingPulsedFreqModFactor();
		this->field->BSDF=dynamic_cast<BDressingPulsedFreqModFactor*>(BDPFMF);		
		this->field->BSDF->fm=parameters.SDparam[3];
  		this->field->BSDF->dw1=parameters.SDparam[4];
  		this->field->BSDF->dw2=parameters.SDparam[5];
  		this->field->BSDF->deltat1=parameters.SDparam[6];
      this->field->BSDF->deltat2=parameters.SDparam[7];
		std::cout<<"Spin Dressing active with Pulsed Frequency Modulation"<<endl;
    std::cout<<this->field->BSDF->deltat2<<endl;
	}
  else if((int)parameters.SDparam[0]==6){
    BDressingFuncModFactor *BDFMF= new BDressingFuncModFactor();  
    this->field->BSDF=dynamic_cast<BDressingFuncModFactor*>(BDFMF);
    this->field->BSDF->fm=parameters.SDparam[3];
    this->field->BSDF->wrf_amp=parameters.SDparam[4];
    this->field->BSDF->scale1=parameters.SDparam[5];
    this->field->BSDF->scale2=parameters.SDparam[6];
    this->field->BSDF->phi_mod=parameters.SDparam[7];
    std::cout<<"Spin Dressing active with Function Frequency Modulation"<<endl;
  }
  else if((int)parameters.SDparam[0]==7){
    BDressingPulsedFreqModFactor2* BDPFMF2 = new BDressingPulsedFreqModFactor2();
    this->field->BSDF=dynamic_cast<BDressingPulsedFreqModFactor2*>(BDPFMF2);    
    this->field->BSDF->fm=parameters.SDparam[3];
      this->field->BSDF->wrf_amp=parameters.SDparam[4];
      this->field->BSDF->scale=parameters.SDparam[5];
      this->field->BSDF->deltat1=parameters.SDparam[6];
      this->field->BSDF->deltat2=parameters.SDparam[7];
    std::cout<<"Spin Dressing active with Pulsed Frequency Modulation 2"<<endl;
    std::cout<<this->field->BSDF->deltat2<<endl;
  }
	else
		for(int i=0;i<100; i++) std::cout<<"The dressing Modulation type is WRONG, Unpredictable behavior imminent"<<endl;
		 		
	
	//parameters all the dressing classes need. 
	this->field->BSDF->wrf=parameters.SDparam[1];
	this->field->BSDF->phi=parameters.SDparam[2];
  this->field->BSDF->noise=parameters.noise;
  this->field->BSDF->pulse=parameters.Pulseparam[0];
  ///Ezra Webb addition of initial pulses!!
  //C Swank addition of noise to pulses! via cmsInterpnoiseGen.C, " ".h
  // noise must be predetermined otherwise no convergence in the RK5 integration! 
  if ((int)parameters.Pulseparam[0]==1){ 
    this->field->BSDF->w_p=parameters.Pulseparam[1];
    this->field->BSDF->Bscale=parameters.Pulseparam[2];
    this->field->BSDF->T_p=parameters.Pulseparam[3];
    this->field->BSDF->phi_p=parameters.Pulseparam[4];
    this->field->BSDF->T_pause=parameters.Pulseparam[5];
    std::cout<<"Initial Pulse with fixed frequency"<<endl;
    
  } 

  else if ((int)parameters.Pulseparam[0]==2){ 
    this->field->BSDF->width=parameters.Pulseparam[1];
    this->field->BSDF->Bscale=parameters.Pulseparam[2];
    this->field->BSDF->T_p=parameters.Pulseparam[3];
    this->field->BSDF->T_crop=parameters.Pulseparam[4];
    this->field->BSDF->T_pause=parameters.Pulseparam[5];
    this->field->BSDF->rscale=parameters.Pulseparam[6];

    std::cout<<"Initial Pulse: sech type"<<endl;
    
  } 

    else if ((int)parameters.Pulseparam[0]==3){ 
    this->field->BSDF->width=parameters.Pulseparam[1];
    this->field->BSDF->Bscale=parameters.Pulseparam[2];
    this->field->BSDF->T_p=parameters.Pulseparam[3];
    this->field->BSDF->T_crop=parameters.Pulseparam[4];
    this->field->BSDF->T_pause=parameters.Pulseparam[5];
    this->field->BSDF->rscale=parameters.Pulseparam[6];
    int tempinterpNum[3]={(int)parameters.Pulseparam[11],(int)parameters.Pulseparam[12],(int)parameters.Pulseparam[13]};
    double tempinterpstart[1]={(double)parameters.Pulseparam[9]};
    double tempinterpstep[1]={(double)parameters.Pulseparam[10]};
    this->field->BSDF->interpnoise=new cmsInterp(tempinterpstart,tempinterpstep,tempinterpNum,(int)parameters.Pulseparam[7],5000,0);
    this->field->BSDF->interpnoise->whitenoiseGen((double)parameters.Pulseparam[8]);
    //double testtesttime=0.0001;
    //this->field->BSDF->interpnoise->testtest=testtesttime;
    //double temptesttime[1]={testtesttime};
    //std::cout<<"What going on??? "<<this->field->BSDF->interpnoise->interp1D(temptesttime)<<"  \n";
    //for(int i=0; i<100; i++) {
      //double temptesttime[1]={(double)i*testtesttime};
      //std::cout<<this->field->BSDF->interpnoise->data[10]<<"  \n";
    //}
    std::cout<<"Initial Pulse: sech type with noise"<<endl;
    }
    else if ((int)parameters.Pulseparam[0]==4){ 
    this->field->BSDF->width=parameters.Pulseparam[1];
    this->field->BSDF->Bscale=parameters.Pulseparam[2];
    this->field->BSDF->T_p=parameters.Pulseparam[3];
    this->field->BSDF->T_crop=parameters.Pulseparam[4];
    this->field->BSDF->T_pause=parameters.Pulseparam[5];
    //this->field->BSDF->rscale=parameters.Pulseparam[6];
    int tempinterpNum[3]={(int)parameters.Pulseparam[11],1,1};
    double tempinterpstart[1]={(double)parameters.Pulseparam[9]};
    double tempinterpstep[1]={(double)parameters.Pulseparam[10]};
    double tempb1interpstart[1]={(double)parameters.Pulseparam[12]};
    double tempb1interpstep[1]={(double)parameters.Pulseparam[13]};
    int tempb1interpNum[3]={(int)parameters.Pulseparam[6],1,1};
    int filenum=(int)parameters.Pulseparam[1];
    //std::string tempb1path ("/data1/cmswank/spin-sim-xliu/BField/B1Pulse.dat");
    //its already hard codded!. 
    this->field->BSDF->interpnoise= new cmsInterp(tempinterpstart,tempinterpstep,tempinterpNum,(int)parameters.Pulseparam[7],5000,0);
    this->field->BSDF->interpnoise->whitenoiseGen((double)parameters.Pulseparam[8]);
    this->field->BSDF->b1Pulse = new cmsB1PulseInterp(tempb1interpstart,tempb1interpstep,tempb1interpNum,(int)filenum);
    //double temptesttime[1]={0.25};
    //std::cout<<"wtf? again? "<<this->field->BSDF->b1Pulse->interp1D(temptesttime)<<"\n";;
    //double testtesttime=0.0001;
    //this->field->BSDF->interpnoise->testtest=testtesttime;
    //
    //std::cout<<"What going on??? "<<this->field->BSDF->interpnoise->interp1D(temptesttime)<<"  \n";
    //for(int i=0; i<100; i++) {
      //double temptesttime[1]={(double)i*testtesttime};
      //std::cout<<this->field->BSDF->interpnoise->data[10]<<"  \n";
    //}
    //int pooperdoop;
    //std::cin>>pooperdoop;
    std::cout<<"Initial Pulse: Numerically tailored with noise"<<endl;
    }

    else if ((int)parameters.Pulseparam[0]==5){ 
    this->field->BSDF->width=parameters.Pulseparam[1];
    this->field->BSDF->Bscale=parameters.Pulseparam[2];
    this->field->BSDF->T_p=parameters.Pulseparam[3];
    this->field->BSDF->T_crop=parameters.Pulseparam[4];
    this->field->BSDF->T_pause=parameters.Pulseparam[5];
    //this->field->BSDF->rscale=parameters.Pulseparam[6];
    int tempinterpNum[3]={(int)parameters.Pulseparam[11],1,1};
    double tempinterpstart[1]={(double)parameters.Pulseparam[9]};
    double tempinterpstep[1]={(double)parameters.Pulseparam[10]};
    double tempb1interpstart[1]={(double)parameters.Pulseparam[12]};
    double tempb1interpstep[1]={(double)parameters.Pulseparam[13]};
    int tempb1interpNum[3]={(int)parameters.Pulseparam[6],1,1};
    int filenum=(int)parameters.Pulseparam[1];
    //std::string tempb1path ("/data1/cmswank/spin-sim-xliu/BField/B1Pulse.dat");
    //its already hard codded!. 
    this->field->BSDF->interpnoise= new cmsInterp(tempinterpstart,tempinterpstep,tempinterpNum,(int)parameters.Pulseparam[7],5000,0);
    this->field->BSDF->interpnoise->whitenoiseGen((double)parameters.Pulseparam[8]);
    this->field->BSDF->b1Pulse = new cmsB1PulseInterp(tempb1interpstart,tempb1interpstep,tempb1interpNum,(int)filenum);
    //double temptesttime[1]={0.};
    
    //std::cout<<"wtf? again? "<<this->field->BSDF->b1Pulse->interp1D(temptesttime)<<"\n";;
    //double testtesttime=0.0001;
    //this->field->BSDF->interpnoise->testtest=testtesttime;
    //
    //std::cout<<"What going on??? "<<this->field->BSDF->interpnoise->interp1D(temptesttime)<<"  \n";
    //for(int i=0; i<100; i++) {
      //double temptesttime[1]={(double)i*testtesttime};
      //std::cout<<this->field->BSDF->interpnoise->data[10]<<"  \n";
    //}
    //int pooperdoop;
    //std::cin>>pooperdoop;
    std::cout<<"Initial Pulse: Robust Dressed Numerically tailored with noise"<<endl;
    }

    else if ((int)parameters.Pulseparam[0]==6){ 
      ///width has been changed to pulse number. 
    this->field->BSDF->width=parameters.Pulseparam[1];
    this->field->BSDF->Bscale=parameters.Pulseparam[2];
    this->field->BSDF->T_p=parameters.Pulseparam[3];
    this->field->BSDF->T_crop=parameters.Pulseparam[4];
    this->field->BSDF->T_pause=parameters.Pulseparam[5];
    
//[0]=6, [1]=width and or filenum, [2]=Bscale, [3]=T_p, [4]= T_crop; [5]=T_pause, [6]= tempb1interpNum. [7]= rng seed, [8]= white noise std. 
    //[9]=tempinterp

    //this->field->BSDF->rscale=parameters.Pulseparam[6];
    int tempinterpNum[3]={(int)parameters.Pulseparam[11],1,1};
    double tempinterpstart[1]={(double)parameters.Pulseparam[9]};
    double tempinterpstep[1]={(double)parameters.Pulseparam[10]};
    double tempb1interpstart[1]={(double)parameters.Pulseparam[12]};
    double tempb1interpstep[1]={(double)parameters.Pulseparam[13]};
    int tempb1interpNum[3]={(int)parameters.Pulseparam[6],1,1};
    int filenum=(int)parameters.Pulseparam[1];
    //std::string tempb1path ("/data1/cmswank/spin-sim-xliu/BField/B1Pulse.dat");
    //its already hard codded!. 
    this->field->BSDF->interpnoise= new cmsInterp(tempinterpstart,tempinterpstep,tempinterpNum,(int)parameters.Pulseparam[7],5000,0);
    this->field->BSDF->interpnoise->whitenoiseGen((double)parameters.Pulseparam[8]);
    this->field->BSDF->b1Pulse = new cmsB1PulseInterp(tempb1interpstart,tempb1interpstep,tempb1interpNum,(int)filenum);
    //double temptesttime[1]={0.25};
    //std::cout<<"wtf? again? "<<this->field->BSDF->b1Pulse->interp1D(temptesttime)<<"\n";;
    //double testtesttime=0.0001;
    //this->field->BSDF->interpnoise->testtest=testtesttime;
    //
    //std::cout<<"What going on??? "<<this->field->BSDF->interpnoise->interp1D(temptesttime)<<"  \n";
    //for(int i=0; i<100; i++) {
      //double temptesttime[1]={(double)i*testtesttime};
      //std::cout<<this->field->BSDF->interpnoise->data[10]<<"  \n";
    //}
    //int pooperdoop;
    //std::cin>>pooperdoop;
    std::cout<<"Pulse is a Dressing pulse found numerically with filtered noise"<<endl;
    }
  else if ((int)parameters.Pulseparam[0]==7){ 
      ///width has been changed to pulse number. 
    this->field->BSDF->width=parameters.Pulseparam[1];
    this->field->BSDF->Bscale=parameters.Pulseparam[2]; //use as Bscale->Wamp [2],// use as T-P->a [3],//use as b T_crop->b [4],//use as T_pause->n [5];
    this->field->BSDF->T_p=parameters.Pulseparam[3];    // use as T-P->a [3]
    this->field->BSDF->T_crop=parameters.Pulseparam[4]; //use as b T_crop->b [4]
    this->field->BSDF->T_pause=parameters.Pulseparam[5]; //use as T_pause->n [5];
    
    //this->field->BSDF->rscale=parameters.Pulseparam[6];
    int tempinterpNum[3]={(int)parameters.Pulseparam[11],1,1};
    double tempinterpstart[1]={(double)parameters.Pulseparam[9]};
    double tempinterpstep[1]={(double)parameters.Pulseparam[10]};
    double highpass={(double)parameters.Pulseparam[13]};
    double cutoff={(double)parameters.Pulseparam[12]};
    int tempb1interpNum[3]={(int)parameters.Pulseparam[6],1,1};
    int filenum=(int)parameters.Pulseparam[1];
    //std::string tempb1path ("/data1/cmswank/spin-sim-xliu/BField/B1Pulse.dat");
    //its already hard codded!. 
    this->field->BSDF->interpnoise= new cmsInterp(tempinterpstart,tempinterpstep,tempinterpNum,(int)parameters.Pulseparam[7],(int)cutoff,(int)highpass);
    this->field->BSDF->interpnoise->whitenoiseGen((double)parameters.Pulseparam[8]);
    //this->field->BSDF->b1Pulse = new cmsB1PulseInterp(tempb1interpstart,tempb1interpstep,tempb1interpNum,(int)filenum);
    //double temptesttime[1]={0.25};
    //std::cout<<"wtf? again? "<<this->field->BSDF->b1Pulse->interp1D(temptesttime)<<"\n";;
    //double testtesttime=0.0001;
    //this->field->BSDF->interpnoise->testtest=testtesttime;
    //
    //std::cout<<"What going on??? "<<this->field->BSDF->interpnoise->interp1D(temptesttime)<<"  \n";
    //for(int i=0; i<100; i++) {
      //double temptesttime[1]={(double)i*testtesttime};
      //std::cout<<this->field->BSDF->interpnoise->data[10]<<"  \n";
    //}
    //int pooperdoop;
    //std::cin>>pooperdoop;
    std::cout<<"Adding Splined Filtered Noise to spin dressing"<<endl;
    }




	if (parameters.field_formula[0])
		this->field->add_field_formulaSD('x', parameters.field_formula[0],this->field->BSDF);
	if (parameters.field_formula[1])
		this->field->add_field_formulaSD('y', parameters.field_formula[1],this->field->BSDF);
	if (parameters.field_formula[2])
		this->field->add_field_formulaSD('z', parameters.field_formula[2],this->field->BSDF);

	}
  else{
  if (parameters.field_formula[0])
    this->field->add_field_formula('x', parameters.field_formula[0]);
  if (parameters.field_formula[1])
    this->field->add_field_formula('y', parameters.field_formula[1]);
  if (parameters.field_formula[2])
    this->field->add_field_formula('z', parameters.field_formula[2]);
	}
  



  // Diffusion
  this->boundary->SetDiffusion(parameters.diffusion);
  // Gravity
  if (parameters.gravity != 0.) {
    this->boundary->SetGravity(parameters.gravity);
  }

  // Initialize factory
  if (!factory) {
    factory = new ParticleFactory();
  }
  factory->parameters = parameters;
  factory->field = field;
  factory->boundary = boundary;
  factory->scattering = scattering;
  for (int i=0; i < 4; i++) {
    factory->vDistribution[i] = vDistribution[i];
  }
}

/**
 * Initialize TFile for storage of data
 */
void Run::InitializeTFile() {
  if (!file) {
    file = new TFile(filename, "RECREATE");

    parameters.SaveParameters();

    dataTree->SetDirectory(file);
    dataTree->Write();
    parameters.tree->Write();

    randomPars->Write("rnd");
    if (parameters.field_formula[0]) {
      parameters.field_formula[0]->Write("field_formula_x");
    }
    if (parameters.field_formula[1]) {
      parameters.field_formula[1]->Write("field_formula_y");
    }
    if (parameters.field_formula[2]) {
      parameters.field_formula[2]->Write("field_formula_z");
    }
  }
}



/**
 * Write data to TFile
 */
void Run::WriteTFile() {
  if (file) {
    file->Write(NULL, TObject::kWriteDelete);
  }
  else {
    InitializeTFile();
  }
}

/**
 * Create data TTree
 */
TTree *Run::createTree(const char* name) {

  TTree *tree = new TTree("dat", "Neutron simulation");
  neutronData.BranchTTree(tree);
  return tree;
}

void Run::SaveNeutronRecord(Neutron *n1, int i) {
  // Save the record in the ntuple (TTree)
  neutronData.bounces[i] = n1->GetBounces();
  neutronData.scat[i] = n1->nScat; // scatterings
  Double_t temp_time = neutronData.time[i] = n1->GetLifetime();
  Double_t *temp = n1->GetPosition();
  neutronData.posx[i] = temp[0];
  neutronData.posy[i] = temp[1];
  neutronData.posz[i] = temp[2];
  temp = n1->GetVelocity();
  neutronData.velx[i] = temp[0];
  neutronData.vely[i] = temp[1];
  neutronData.velz[i] = temp[2];
  neutronData.speed[i] = Vector::Norm(temp);
  temp = n1->GetSpin();
  neutronData.spinx[i] = temp[0];
  neutronData.spiny[i] = temp[1];
  neutronData.spinz[i] = temp[2];

  BFieldVars temp_vars(temp_time, n1->GetPosition(), n1->GetVelocity());
  Double_t temp4[3];
  field->getField(temp4, temp_vars);
  //if (i==3)
    //std::cout<<temp4[2]<<std::endl;
    
  neutronData.fieldx[i] = temp4[0];
  neutronData.fieldy[i] = temp4[1];
  neutronData.fieldz[i] = temp4[2];

  neutronData.steps[i] = n1->GetProp();
  neutronData.backsteps[i] = n1->GetBackSteps();
};

/**
 * Run a single particle (Neutron)
 */
Int_t Run::runNeutron() {

  Neutron *n1 = factory->newParticle();

  //std::cout<<"whats the field vector size? "<<factory->field->fields.size()<<std::endl;
  //double testpos[3]={0.,0.,0.};
  //double testvel[3]={0.,0.,0.,};
  //BFieldVars testvars(1./3845./2.,testpos,testvel);
  //double *Btest=new double[3];
  
  //factory->field->getField(Btest,testvars);
  
  //std::cout<<"What is about the factory field value? "<<Btest[2]<<std::endl;
  //factory->field->fields.erase(factory->field->fields.begin(),factory->field->fields.begin()+factory->field->fields.size());
 // factory->field->getField(Btest,testvars);
  //std::cout<<"how about now? "<<Btest[2]<<std::endl;
  //Double_t gamE=n1->GetGammaE();
  //Double_t nEDM=n1->edm;
  //Double_t gam=n1->GetGamma();
  //std::cout<<gamE<<" ";
  //std::cout<<nEDM<<" ";

  Double_t nextTime;
  Int_t nBounces = 0;
  Int_t collision = 0;

  Int_t i; // equally spaced samples to take
  // in addition to t=0 which is ALWAYS stored at the zeroth place
  // spacing is done so that the last time "totalTime" is the last entry

  Double_t *temp = 0;

  neutronData.n = parameters.nBins + 1;
  neutronData.id = nCurrent++;

  static Double_t *randBinTimes = 0;
  static Int_t *randBinOrder = 0;
  if (parameters.randBins > 0 && !randBinTimes) {
    // the times are randomly distributed over totalTime
    // NB the bins are sorted by time
    // so the random time distribution in each bin is not uniform 
    randBinTimes = new Double_t[parameters.nBins];
    randBinOrder = new Int_t[parameters.nBins];
    for (i = 0; i < parameters.nBins; i++) {
      randBinTimes[i] = gRandom->Rndm() * (parameters.totalTime - parameters.offsetTime) + parameters.offsetTime;
    }
    TMath::Sort(parameters.nBins, randBinTimes, randBinOrder, 0);
    // descending order sort
  }

  nBounces = 0;nextTime = 0.;
  for (i = 0; i <= parameters.nBins; i++) {
    if (i > 0 && n1->GetLifetime() < nextTime)
      cout << "advanced not hard enough i=" << i << endl;
    
    if (parameters.randBins > 0 && i < parameters.nBins) {
      nextTime = randBinTimes[randBinOrder[i]];
    }
    else {
      // calculate the next snapshot time
      nextTime = parameters.offsetTime + ((Double_t)i + 1.) * (parameters.totalTime - parameters.offsetTime) / ((Double_t)parameters.nBins);
    }

    SaveNeutronRecord(n1, i);

    // Save the first position as "collision"
    if (i == 0) {
      temp = n1->GetPosition();
      neutronData.bx[0] = temp[0];
      neutronData.by[0] = temp[1];
      neutronData.bz[0] = temp[2];
      neutronData.bkind[0] = 0;
    }

    // When done saving, exit the for loop
    if (i == parameters.nBins)
      break;

    AdvanceOption advance_mode;
    if (nBounces < 200) {
      advance_mode = ADVANCE_STOP_ON_BOUNCE;
    }
    else {
      advance_mode = ADVANCE_NO_STOP;
    }

    // Propagate the particle
    while (n1->GetLifetime() < nextTime) {
      Double_t delta_t = nextTime - n1->GetLifetime();
      Double_t advance_ret = n1->Advance(delta_t, advance_mode);
      //gam=n1->GetGamma();
      
      //if (gam<-193247173) std::cout<<gam<<"  ";
      

      if (advance_ret < 0) {
	cout << "Error: Advance returned time < 0" << endl;
	delete n1;
	return -1;
      }
      // In case of a collision
      if (collision > 0) {
	nBounces++;
	if (nBounces < 200) {
	  temp = n1->GetPosition();
	  neutronData.bx[nBounces] = temp[0];
	  neutronData.by[nBounces] = temp[1];
	  neutronData.bz[nBounces] = temp[2];
	  neutronData.bkind[nBounces] = collision;
	} 
      }
      
    } // end while loop    
      
  } // end for loop

  // This is where it fills the data tree. We need to remake the data
  dataTree->Fill();

  if (parameters.randBins == 2) {
    // set randBins = 2; to randomize the bins for each particle
    delete[] randBinTimes;
    delete[] randBinOrder;
    randBinTimes = 0;
  }

  delete n1;

  return 0;
}





Int_t Run::runNeutronQuiet() {

  Neutron *n1 = factory->newParticle();
  //std::cout<<"whats the field vector size? "<<factory->field->fields.size()<<std::endl;
  //double testpos[3]={0.,0.,0.};
  //double testvel[3]={0.,0.,0.,};
  //BFieldVars testvars(1./3845./2.,testpos,testvel);
  //double *Btest=new double[3];
  
  //factory->field->getField(Btest,testvars);
  
  //std::cout<<"What is about the factory field value? "<<Btest[2]<<std::endl;
  //factory->field->fields.erase(factory->field->fields.begin(),factory->field->fields.begin()+factory->field->fields.size());
 // factory->field->getField(Btest,testvars);
  //std::cout<<"how about now? "<<Btest[2]<<std::endl;
  
  //Double_t gam=n1->GetGamma();
  //std::cout<<gam<<" ";

  Double_t nextTime;
  Int_t nBounces = 0;
  Int_t collision = 0;

  Int_t i; // equally spaced samples to take
  // in addition to t=0 which is ALWAYS stored at the zeroth place
  // spacing is done so that the last time "totalTime" is the last entry

  Double_t *temp = 0;

  neutronData.n = parameters.nBins + 1;
  neutronData.id = nCurrent++;

  static Double_t *randBinTimes = 0;
  static Int_t *randBinOrder = 0;
  if (parameters.randBins > 0 && !randBinTimes) {
    // the times are randomly distributed over totalTime
    // NB the bins are sorted by time
    // so the random time distribution in each bin is not uniform 
    randBinTimes = new Double_t[parameters.nBins];
    randBinOrder = new Int_t[parameters.nBins];
    for (i = 0; i < parameters.nBins; i++) {
      randBinTimes[i] = gRandom->Rndm() * (parameters.totalTime - parameters.offsetTime) + parameters.offsetTime;
    }
    TMath::Sort(parameters.nBins, randBinTimes, randBinOrder, 0);
    // descending order sort
  }

  nBounces = 0;nextTime = 0.;
  for (i = 0; i <= parameters.nBins; i++) {
    if (i > 0 && n1->GetLifetime() < nextTime)
      cout << "advanced not hard enough i=" << i << endl;
    
    if (parameters.randBins > 0 && i < parameters.nBins) {
      nextTime = randBinTimes[randBinOrder[i]];
    }
    else {
      // calculate the next snapshot time
      nextTime = parameters.offsetTime + ((Double_t)i + 1.) * (parameters.totalTime - parameters.offsetTime) / ((Double_t)parameters.nBins);
    }

    SaveNeutronRecord(n1, i);

    // Save the first position as "collision"
    if (i == 0) {
      temp = n1->GetPosition();
      neutronData.bx[0] = temp[0];
      neutronData.by[0] = temp[1];
      neutronData.bz[0] = temp[2];
      neutronData.bkind[0] = 0;
    }

    // When done saving, exit the for loop
    if (i == parameters.nBins)
      break;

    AdvanceOption advance_mode;
    if (nBounces < 200) {
      advance_mode = ADVANCE_STOP_ON_BOUNCE;
    }
    else {
      advance_mode = ADVANCE_NO_STOP;
    }

    // Propagate the particle
    while (n1->GetLifetime() < nextTime) {
      Double_t delta_t = nextTime - n1->GetLifetime();
      Double_t advance_ret = n1->Advance(delta_t, advance_mode);
      //gam=n1->GetGamma();
      
      //if (gam<-193247173) std::cout<<gam<<"  ";
      

      if (advance_ret < 0) {
  cout << "Error: Advance returned time < 0" << endl;
  delete n1;
  return -1;
      }
      // In case of a collision
      if (collision > 0) {
  nBounces++;
  if (nBounces < 200) {
    temp = n1->GetPosition();
    neutronData.bx[nBounces] = temp[0];
    neutronData.by[nBounces] = temp[1];
    neutronData.bz[nBounces] = temp[2];
    neutronData.bkind[nBounces] = collision;
  } 
      }
      
    } // end while loop    
      
  } // end for loop

  // This is where it fills the data tree. We need to remake the data
  //dataTree->Fill();

  if (parameters.randBins == 2) {
    // set randBins = 2; to randomize the bins for each particle
    delete[] randBinTimes;
    delete[] randBinOrder;
    randBinTimes = 0;
  }

  delete n1;

  return 0;
}













/**
 * Run all particles in loop
 */
int Run::runAll() {
  int count = 0;
  for (int i=nCurrent; i < parameters.nEntries && !RUN_EXIT_ASAP; i++) {
    if (runNeutron() == 0) {
      cout << "." << flush;
      count++;
    }
    else {
      i--;
      cout << "Discarded one neutron" << endl;
    }
  }
  return count;
}

/**
 * Use the Maxwell-Boltzmann distribution
 * to find speed of atoms
 */
void Run::useBoltzmannDistribution(Double_t temperature) {

  TF1 *f1;

  const Double_t boltzmann = 1.38066e-23; // J K^-1
  const Double_t MHe3 = 2.2 * 3.016029 * 1.66054e-27; // kg !!effective Mass (M*=2.2M)

  cout << "Using Boltzmann (T = " << temperature << " K)" << endl;

  Double_t temp = MHe3 / (2. * temperature * boltzmann);

  f1 = new TF1("boltzmann", "(x^2) * exp(-(x^2) * [0])", 0, 
	       10. * sqrt(1. / temp));
  // limit set to 10 v_p (the most probable velocity)
  // 4.2e-43 integrated probability beyond that velocity

  f1->SetNpx(2000);
  
  f1->SetParameter(0,temp);

  this->vDistribution[0] = f1;


  temp = sqrt(8. * boltzmann * temperature / (MHe3 * 3.14159));

  this->mfp = 3. * (1.6e-4 / pow(temperature, 7)) / temp;

  temp = sqrt(3. * boltzmann * temperature / MHe3);

  this->tau_C = 3. * (1.6e-4 / pow(temperature, 7)) / (temp*temp);

  cout << "MFP = " << this->mfp << " m" << endl;

  this->scattering = new Scattering(1);
  this->scattering->Set(this->tau_C);

  cout << "Scattering class enabled with \\tau_C = " << this->tau_C << endl;

  return;


}

/**
 * Check if the data file already exists
 */
int Run::FileExists(const char *filename) {

  int fileExists = 0;

  ifstream *tempFile = new ifstream();

  tempFile->open(filename);
  if (tempFile->is_open())
    fileExists = 1;
  tempFile->close();

  return fileExists;    

}

/**
 * Parse options on argument line
 */
void Run::parse_argument_line(int argc, char **argv, int *options, char **filename) {

  // Get other arguments  
  int  argi = 2;
  while (argi < argc) {  // Get the other parameters
    // -p parameter_file
    if (argc - argi >= 2 && strcmp(argv[argi], "-p") == 0) {      
      options[0] = 1;
      *filename = argv[argi + 1];
      argi += 2;      
    }
    // -c 
    else if (argc - argi >= 1 && strcmp(argv[argi], "-c") == 0) {
      options[1] = 1;
      argi++;
    }
    // -n #entries
    else if (argc - argi >= 2 && strcmp(argv[argi], "-n") == 0) {
      options[2] = atoi(argv[argi + 1]);
      argi += 2;
    }
    // -novc   (no version check)
    else if (argc - argi >= 1 && strcmp(argv[argi], "--novc") == 0) {
      options[3] = 1;
      argi++;
    }
    // -numerical
    else if (argc - argi >= 2 && strcmp(argv[argi], "--numerical") == 0) {
      options[4] = atoi(argv[argi + 1]);
      argi += 2;
    }
    else {
      cout << "Argument line not understood" << endl;
      exit(0);
    }

  } // end while loop
  

}

int RUN_NEEDS_SAVE = 0;
int RUN_EXIT_ASAP = 0;

/**
 * Signal handler, will quit gracefully if interrupted
 */
void signal_handle_int(int signal) {
  cout << endl;

  if (RUN_EXIT_ASAP) {
    cout << "OK exiting now!" << endl;
    exit(1);
  }

  RUN_EXIT_ASAP = 1;
  if (RUN_NEEDS_SAVE) {
    cout << "Letting file finish saving..." << endl;
  }
  return;
}

/**
 * Set signal handlers
 */
void set_signal_handles() {
  signal(SIGINT, signal_handle_int);
  signal(SIGTERM, signal_handle_int);
}


//CS addition Erase field formula
void Run::EraseFieldFormula(){
  for (int findi=0; findi<3; findi++)
    this->parameters.field_formula[findi]=0;
}

//CS Addition ReloadParameters. 

void Run::ReloadParameters() {

  //neutronData.Initialize(parameters.nBins + 1);

  //if (!dataTree) {
  //  dataTree = createTree(name);
  //}
  //std::cout<<"Reloding Parameters"<<std::endl;

  Double_t gamma = 0;
  switch (parameters.simType) {
  case 0:
    //cout << "Neutron Simulation" << endl;
    gamma = NEUTRON_GAMMA;
    //std::cout<<"neutron gamma: "<<gamma<<std::endl;
    break;
  case 1:
    // cout << "He3 Simulation" << endl;
    gamma = HE3_GAMMA;
    //std::cout<<"3He gamma: "<<gamma<<std::endl;
    useQuiteBoltzmannDistribution(parameters.temperature);
    break;
  }
  if (parameters.simGamma != 0) {
    gamma = parameters.simGamma;
    //cout << "User defined gamma = " << gamma << endl;
  }
  //cout.precision(10);
  //cout << "gamma = " << gamma << endl;

  // make sure that offset time is positive and < totalTime
  if (parameters.offsetTime < 0. || parameters.offsetTime >= parameters.totalTime) {
    //cout << "# invalid offset time provided" << endl;
    parameters.offsetTime = 0.;
  }

  // override the standard gRandom global variable
  // using TRandom3 which is a better generator
  // the generator is seeded with 0 so that it chooses its's 
  // seed based on time??
  
  //cs looking fore seg fault 
  //how can this work?
  //delete gRandom;
  //if (this->randomPars) {
   // gRandom = randomPars;
    //cout << "rnd restored" << endl;
  //}
  //else {
    //if (parameters.seed == 0) {
      //time_t curtime;
      //time(&curtime);
      //parameters.seed = (Int_t)curtime;
    //}
    gRandom = new TRandom3(parameters.seed);
    parameters.seed = gRandom->GetSeed();
    //cout << "seed : " << parameters.seed << endl;
    randomPars = gRandom;
  //}


  // Set Geometry
  /*
  if(boundary) delete boundary;

  switch (parameters.geometry) {
  case 0:
    boundary = new CircleBoundary(parameters.radius); // 2D circle
    break;
  case 1:
    // box
    boundary = new BoxBoundary(parameters.boxLow, parameters.boxHigh);
    break;
  case 2:
    boundary = new SphereBoundary(parameters.radius); // sphere
    break;
  case 3:
    boundary = new Box2DBoundary(parameters.boxLow, parameters.boxHigh);
    break;
  case 4:
    boundary = new CylindricalBoundary(parameters.radius, parameters.boxLow[2], parameters.boxHigh[2]);
    break;
  case 5:
    boundary = new Boundary1D(parameters.boxLow[0], parameters.boxHigh[0]);
    break;
  default:
    boundary = 0;
  }*/

    //remove previous fields, this turned out to be tricky because its a mass of spagetti and a lot are protected. 
  
  //this->field->BSDF->interpnoise->~cmsInterp();

  //this->field->BSDF->b1Pulse->~cmsB1PulseInterp();
  //delete this->field;

  // Set Uniform Fields
    
  this->field = new ParticleField(parameters.B0, parameters.E0);
  if (parameters.uniform_gradient)
    this->field->SetUniformGrad(parameters.uniform_gradient);
  else {
    if (parameters.simple_gradient) {
      field->SetSimpleGrad(parameters.simple_gradient);
    }
    if (parameters.simple_quad_gradient) {
      field->SetSimpleQuadGrad(parameters.simple_quad_gradient);
    }
  }
  //std::cout<<"past simple fiels"<<std::endl;

  // Set additional fields
// Set Interpolation fields (generated by COMSOL or other FEM)
///added by C.Swank, TODO include time dependence
if (parameters.Interpparam[0]){
  
  //time independent? 
  if ((int)parameters.Interpparam[0]==1){ 
    
  double star[3]={(double)parameters.Interpparam[1],(double)parameters.Interpparam[2],(double)parameters.Interpparam[3]};
  double ste[3]={(double)parameters.Interpparam[4],(double)parameters.Interpparam[5],(double)parameters.Interpparam[6]};
  int Ninterp[3]={(int)parameters.Interpparam[7],(int)parameters.Interpparam[8],(int)parameters.Interpparam[9]};
  double nuttin[1]={0.};
    
    //std::cout<<"B Field Interpolation"<<endl;
   //constructor for reference. 
   //BFieldInterp::BFieldInterp(bool fromfile, int dir,double *star, double *ste, int N[3], double *datainput)
   //std::cout<<parameters.InterpDir<<"\n";
    
  if (strcmp(parameters.InterpDir.c_str(),"x")==0)
    {
      this->field->BFIx=new BFieldInterp(true,0,star,ste,Ninterp,nuttin);  
        this->field->add_field_interp('x', this->field->BFIx);
      //std::cout<<"adding Bx interpolation field\n";
      //double tempos[3]={0.3280196619,0.1629835684,-0.8445006203};
      //std::cout<<"test fied "<<this->field->BFIx->interp3D(tempos)<<"\n";
    }
  if (strcmp(parameters.InterpDir.c_str(),"y")==0)
    {
    this->field->BFIy=new BFieldInterp(true,1,star,ste,Ninterp,nuttin);  
    this->field->add_field_interp('y', this->field->BFIy);
    //std::cout<<"adding By interpolation field\n";
    }
  if (strcmp(parameters.InterpDir.c_str(),"z")==0)
  {
    this->field->BFIz=new BFieldInterp(true,2,star,ste,Ninterp,nuttin); 
    this->field->add_field_interp('z', this->field->BFIz);
    //std::cout<<"adding Bz interpolation field\n";
  }
  } 
}


  //std::cout<<"past interp field"<<std::endl;


  if (parameters.SDparam[0]){
  
  if ((int)parameters.SDparam[0]==1){ 
    this->field->BSDF=new BDressingFactor();  
    //std::cout<<"Spin Dressing active without Modulation"<<endl;
    
  } 
  else if((int)parameters.SDparam[0]==2){
    BDressingCosModFactor *BDCMF= new BDressingCosModFactor();  
    this->field->BSDF=dynamic_cast<BDressingCosModFactor*>(BDCMF);
    this->field->BSDF->fm=parameters.SDparam[3];
    this->field->BSDF->wrf_amp=parameters.SDparam[4];
    this->field->BSDF->wrf=parameters.SDparam[1];
    this->field->BSDF->t1=parameters.SDparam[5];
    this->field->BSDF->t2=parameters.SDparam[6];
    this->field->BSDF->nmod=parameters.SDparam[7];
    this->field->BSDF->findFM();
    //std::cout<<"Spin Dressing active with Cosine Frequency Modulation"<<endl;
    
  }
  else if((int)parameters.SDparam[0]==3){
    BDressingCosBModFactor *BDCBMF= new BDressingCosBModFactor(); 
    this->field->BSDF=dynamic_cast<BDressingCosBModFactor*>(BDCBMF);
    this->field->BSDF->fm=parameters.SDparam[3];
    this->field->BSDF->amp=parameters.SDparam[4];
      //std::cout<<"Spin Dressing active with Cosine Field Magnitude Modulation"<<endl;
  
  }
  else if((int)parameters.SDparam[0]==4){
    BDressingPulsedBModFactor* BDPBMF = new BDressingPulsedBModFactor();
    this->field->BSDF=dynamic_cast<BDressingPulsedBModFactor*>(BDPBMF);
    this->field->BSDF->fm=parameters.SDparam[3];
      this->field->BSDF->scale1=parameters.SDparam[4];
      this->field->BSDF->scale2=parameters.SDparam[5];
      this->field->BSDF->deltat=parameters.SDparam[6];
    //std::cout<<"Spin Dressing active with Pulsed Field Magnitude Modulation"<<endl;
  }
  else if((int)parameters.SDparam[0]==5){
    BDressingPulsedFreqModFactor* BDPFMF = new BDressingPulsedFreqModFactor();
    this->field->BSDF=dynamic_cast<BDressingPulsedFreqModFactor*>(BDPFMF);    
    this->field->BSDF->fm=parameters.SDparam[3];
      this->field->BSDF->dw1=parameters.SDparam[4];
      this->field->BSDF->dw2=parameters.SDparam[5];
      this->field->BSDF->deltat1=parameters.SDparam[6];
      this->field->BSDF->deltat2=parameters.SDparam[7];
    //std::cout<<"Spin Dressing active with Pulsed Frequency Modulation"<<endl;
    //std::cout<<this->field->BSDF->deltat2<<endl;
  }
  else if((int)parameters.SDparam[0]==6){
    BDressingFuncModFactor *BDFMF= new BDressingFuncModFactor();  
    this->field->BSDF=dynamic_cast<BDressingFuncModFactor*>(BDFMF);
    this->field->BSDF->fm=parameters.SDparam[3];
    this->field->BSDF->wrf_amp=parameters.SDparam[4];
    this->field->BSDF->scale1=parameters.SDparam[5];
    this->field->BSDF->scale2=parameters.SDparam[6];
    this->field->BSDF->phi_mod=parameters.SDparam[7];
    //std::cout<<"Spin Dressing active with Function Frequency Modulation"<<endl;
  }
  else if((int)parameters.SDparam[0]==7){
    BDressingPulsedFreqModFactor2* BDPFMF2 = new BDressingPulsedFreqModFactor2();
    this->field->BSDF=dynamic_cast<BDressingPulsedFreqModFactor2*>(BDPFMF2);    
    this->field->BSDF->fm=parameters.SDparam[3];
      this->field->BSDF->wrf_amp=parameters.SDparam[4];
      this->field->BSDF->scale=parameters.SDparam[5];
      this->field->BSDF->deltat1=parameters.SDparam[6];
      this->field->BSDF->deltat2=parameters.SDparam[7];
    //std::cout<<"Spin Dressing active with Pulsed Frequency Modulation 2"<<endl;
    //std::cout<<this->field->BSDF->deltat2<<endl;
  }
  else
    for(int i=0;i<100; i++) std::cout<<"The dressing Modulation type is WRONG, Unpredictable behavior imminent"<<endl;
        
  
  //parameters all the dressing classes need. 
  this->field->BSDF->wrf=parameters.SDparam[1];
  this->field->BSDF->phi=parameters.SDparam[2];
  this->field->BSDF->noise=parameters.noise;
  this->field->BSDF->pulse=parameters.Pulseparam[0];
  

  //std::cout<<"past the SD field"<<std::endl;

  ///Ezra Webb addition of initial pulses!!
  //C Swank addition of noise to pulses! via cmsInterpnoiseGen.C, " ".h
  // noise must be predetermined otherwise no convergence in the RK5 integration! 
  if ((int)parameters.Pulseparam[0]==1){ 
    this->field->BSDF->w_p=parameters.Pulseparam[1];
    this->field->BSDF->Bscale=parameters.Pulseparam[2];
    this->field->BSDF->T_p=parameters.Pulseparam[3];
    this->field->BSDF->phi_p=parameters.Pulseparam[4];
    this->field->BSDF->T_pause=parameters.Pulseparam[5];
    //std::cout<<"Initial Pulse with fixed frequency"<<endl;
    
  } 

  else if ((int)parameters.Pulseparam[0]==2){ 
    this->field->BSDF->width=parameters.Pulseparam[1];
    this->field->BSDF->Bscale=parameters.Pulseparam[2];
    this->field->BSDF->T_p=parameters.Pulseparam[3];
    this->field->BSDF->T_crop=parameters.Pulseparam[4];
    this->field->BSDF->T_pause=parameters.Pulseparam[5];
    this->field->BSDF->rscale=parameters.Pulseparam[6];

    //std::cout<<"Initial Pulse: sech type"<<endl;
    
  } 

    else if ((int)parameters.Pulseparam[0]==3){ 
    this->field->BSDF->width=parameters.Pulseparam[1];
    this->field->BSDF->Bscale=parameters.Pulseparam[2];
    this->field->BSDF->T_p=parameters.Pulseparam[3];
    this->field->BSDF->T_crop=parameters.Pulseparam[4];
    this->field->BSDF->T_pause=parameters.Pulseparam[5];
    this->field->BSDF->rscale=parameters.Pulseparam[6];
    int tempinterpNum[3]={(int)parameters.Pulseparam[11],(int)parameters.Pulseparam[12],(int)parameters.Pulseparam[13]};
    double tempinterpstart[1]={(double)parameters.Pulseparam[9]};
    double tempinterpstep[1]={(double)parameters.Pulseparam[10]};
    this->field->BSDF->interpnoise=new cmsInterp(tempinterpstart,tempinterpstep,tempinterpNum,(int)parameters.Pulseparam[7],5000,0);
    this->field->BSDF->interpnoise->whitenoiseGen((double)parameters.Pulseparam[8]);
    //double testtesttime=0.0001;
    //this->field->BSDF->interpnoise->testtest=testtesttime;
    //double temptesttime[1]={testtesttime};
    //std::cout<<"What going on??? "<<this->field->BSDF->interpnoise->interp1D(temptesttime)<<"  \n";
    //for(int i=0; i<100; i++) {
      //double temptesttime[1]={(double)i*testtesttime};
      //std::cout<<this->field->BSDF->interpnoise->data[10]<<"  \n";
    //}
    //std::cout<<"Initial Pulse: sech type with noise"<<endl;
    }
    else if ((int)parameters.Pulseparam[0]==4){ 
    this->field->BSDF->width=parameters.Pulseparam[1];
    this->field->BSDF->Bscale=parameters.Pulseparam[2];
    this->field->BSDF->T_p=parameters.Pulseparam[3];
    this->field->BSDF->T_crop=parameters.Pulseparam[4];
    this->field->BSDF->T_pause=parameters.Pulseparam[5];
    //this->field->BSDF->rscale=parameters.Pulseparam[6];
    int tempinterpNum[3]={(int)parameters.Pulseparam[11],1,1};
    double tempinterpstart[1]={(double)parameters.Pulseparam[9]};
    double tempinterpstep[1]={(double)parameters.Pulseparam[10]};
    double tempb1interpstart[1]={(double)parameters.Pulseparam[12]};
    double tempb1interpstep[1]={(double)parameters.Pulseparam[13]};
    int tempb1interpNum[3]={(int)parameters.Pulseparam[6],1,1};
    int filenum=(int)parameters.Pulseparam[1];
    //std::string tempb1path ("/data1/cmswank/spin-sim-xliu/BField/B1Pulse.dat");
    //its already hard codded!. 
    this->field->BSDF->interpnoise= new cmsInterp(tempinterpstart,tempinterpstep,tempinterpNum,(int)parameters.Pulseparam[7],5000,0);
    this->field->BSDF->interpnoise->whitenoiseGen((double)parameters.Pulseparam[8]);
    this->field->BSDF->b1Pulse = new cmsB1PulseInterp(tempb1interpstart,tempb1interpstep,tempb1interpNum,(int)filenum);
    //double temptesttime[1]={0.25};
    //std::cout<<"wtf? again? "<<this->field->BSDF->b1Pulse->interp1D(temptesttime)<<"\n";;
    //double testtesttime=0.0001;
    //this->field->BSDF->interpnoise->testtest=testtesttime;
    //
    //std::cout<<"What going on??? "<<this->field->BSDF->interpnoise->interp1D(temptesttime)<<"  \n";
    //for(int i=0; i<100; i++) {
      //double temptesttime[1]={(double)i*testtesttime};
      //std::cout<<this->field->BSDF->interpnoise->data[10]<<"  \n";
    //}
    //int pooperdoop;
    //std::cin>>pooperdoop;
    //std::cout<<"Initial Pulse: Numerically tailored with noise"<<endl;
    }

    else if ((int)parameters.Pulseparam[0]==5){ 
    this->field->BSDF->width=parameters.Pulseparam[1];
    this->field->BSDF->Bscale=parameters.Pulseparam[2];
    this->field->BSDF->T_p=parameters.Pulseparam[3];
    this->field->BSDF->T_crop=parameters.Pulseparam[4];
    this->field->BSDF->T_pause=parameters.Pulseparam[5];
    //this->field->BSDF->rscale=parameters.Pulseparam[6];
    int tempinterpNum[3]={(int)parameters.Pulseparam[11],1,1};
    double tempinterpstart[1]={(double)parameters.Pulseparam[9]};
    double tempinterpstep[1]={(double)parameters.Pulseparam[10]};
    double tempb1interpstart[1]={(double)parameters.Pulseparam[12]};
    double tempb1interpstep[1]={(double)parameters.Pulseparam[13]};
    int tempb1interpNum[3]={(int)parameters.Pulseparam[6],1,1};
    int filenum=(int)parameters.Pulseparam[1];
    //std::string tempb1path ("/data1/cmswank/spin-sim-xliu/BField/B1Pulse.dat");
    //its already hard codded!. 
    this->field->BSDF->interpnoise= new cmsInterp(tempinterpstart,tempinterpstep,tempinterpNum,(int)parameters.Pulseparam[7],5000,0);
    this->field->BSDF->interpnoise->whitenoiseGen((double)parameters.Pulseparam[8]);
    this->field->BSDF->b1Pulse = new cmsB1PulseInterp(tempb1interpstart,tempb1interpstep,tempb1interpNum,(int)filenum);
    //double temptesttime[1]={0.};
    
    //std::cout<<"wtf? again? "<<this->field->BSDF->b1Pulse->interp1D(temptesttime)<<"\n";;
    //double testtesttime=0.0001;
    //this->field->BSDF->interpnoise->testtest=testtesttime;
    //
    //std::cout<<"What going on??? "<<this->field->BSDF->interpnoise->interp1D(temptesttime)<<"  \n";
    //for(int i=0; i<100; i++) {
      //double temptesttime[1]={(double)i*testtesttime};
      //std::cout<<this->field->BSDF->interpnoise->data[10]<<"  \n";
    //}
    //int pooperdoop;
    //std::cin>>pooperdoop;
    //std::cout<<"Initial Pulse: Robust Dressed Numerically tailored with noise"<<endl;
    }

    else if ((int)parameters.Pulseparam[0]==6){ 
      ///width has been changed to pulse number. 
    this->field->BSDF->width=parameters.Pulseparam[1];
    this->field->BSDF->Bscale=parameters.Pulseparam[2];
    this->field->BSDF->T_p=parameters.Pulseparam[3];
    this->field->BSDF->T_crop=parameters.Pulseparam[4];
    this->field->BSDF->T_pause=parameters.Pulseparam[5];
    
    //this->field->BSDF->rscale=parameters.Pulseparam[6];
    int tempinterpNum[3]={(int)parameters.Pulseparam[11],1,1};
    double tempinterpstart[1]={(double)parameters.Pulseparam[9]};
    double tempinterpstep[1]={(double)parameters.Pulseparam[10]};
    double tempb1interpstart[1]={(double)parameters.Pulseparam[12]};
    double tempb1interpstep[1]={(double)parameters.Pulseparam[13]};
    int tempb1interpNum[3]={(int)parameters.Pulseparam[6],1,1};
    int filenum=(int)parameters.Pulseparam[1];
    //std::string tempb1path ("/data1/cmswank/spin-sim-xliu/BField/B1Pulse.dat");
    //its already hard codded!. 
    this->field->BSDF->interpnoise= new cmsInterp(tempinterpstart,tempinterpstep,tempinterpNum,(int)parameters.Pulseparam[7],5000,0);
    this->field->BSDF->interpnoise->whitenoiseGen((double)parameters.Pulseparam[8]);
    this->field->BSDF->b1Pulse = new cmsB1PulseInterp(tempb1interpstart,tempb1interpstep,tempb1interpNum,(int)filenum);
    //std::cout<<this->field->BSDF->b1Pulse.filenum<<std::endl;
    //double temptesttime[1]={0.25};
    //std::cout<<"wtf? again? "<<this->field->BSDF->b1Pulse->interp1D(temptesttime)<<"\n";;
    //double testtesttime=0.0001;
    //this->field->BSDF->interpnoise->testtest=testtesttime;
    //
    //std::cout<<"What going on??? "<<this->field->BSDF->interpnoise->interp1D(temptesttime)<<"  \n";
    //for(int i=0; i<100; i++) {
      //double temptesttime[1]={(double)i*testtesttime};
      //std::cout<<this->field->BSDF->interpnoise->data[10]<<"  \n";
    //}
    //int pooperdoop;
    //std::cin>>pooperdoop;
    //std::cout<<"Pulse is a Dressing pulse found numerically with filtered noise"<<endl;
    }
    else if ((int)parameters.Pulseparam[0]==7){ 
      ///width has been changed to pulse number. 
    this->field->BSDF->width=parameters.Pulseparam[1];
    this->field->BSDF->Bscale=parameters.Pulseparam[2];
    this->field->BSDF->T_p=parameters.Pulseparam[3];
    this->field->BSDF->T_crop=parameters.Pulseparam[4];
    this->field->BSDF->T_pause=parameters.Pulseparam[5];
    
    //this->field->BSDF->rscale=parameters.Pulseparam[6];
    int tempinterpNum[3]={(int)parameters.Pulseparam[11],1,1};
    double tempinterpstart[1]={(double)parameters.Pulseparam[9]};
    double tempinterpstep[1]={(double)parameters.Pulseparam[10]};
    double highpass={(double)parameters.Pulseparam[13]};
    double cutoff={(double)parameters.Pulseparam[12]};
    int tempb1interpNum[3]={(int)parameters.Pulseparam[6],1,1};
    int filenum=(int)parameters.Pulseparam[1];
    //std::string tempb1path ("/data1/cmswank/spin-sim-xliu/BField/B1Pulse.dat");
    //its already hard codded!. 
    this->field->BSDF->interpnoise= new cmsInterp(tempinterpstart,tempinterpstep,tempinterpNum,(int)parameters.Pulseparam[7],(int)cutoff,(int)highpass);
    //std::cout<<this->field->BSDF->interpnoise->rnd_seed<<std::endl;
    //std::cout<<"HI"<<std::flush;
    this->field->BSDF->interpnoise->whitenoiseGen((double)parameters.Pulseparam[8]);
    //std::cout<<"made it?"<<std::flush;
    //this->field->BSDF->b1Pulse = new cmsB1PulseInterp(tempb1interpstart,tempb1interpstep,tempb1interpNum,(int)filenum);
    //double temptesttime[1]={0.25};
    //std::cout<<"wtf? again? "<<this->field->BSDF->b1Pulse->interp1D(temptesttime)<<"\n";;
    //double testtesttime=0.0001;
    //this->field->BSDF->interpnoise->testtest=testtesttime;
    //
    //std::cout<<"What going on??? "<<this->field->BSDF->interpnoise->interp1D(temptesttime)<<"  \n";
    //for(int i=0; i<100; i++) {
      //double temptesttime[1]={(double)i*testtesttime};
      //std::cout<<this->field->BSDF->interpnoise->data[10]<<"  \n";
    //}
    //int pooperdoop;
    //std::cin>>pooperdoop;
    //std::cout<<"Adding Splined Filtered Noise to spin dressing"<<endl;
    }





  if (parameters.field_formula[0])
    this->field->add_field_formulaSD('x', parameters.field_formula[0],this->field->BSDF);
  if (parameters.field_formula[1])
    this->field->add_field_formulaSD('y', parameters.field_formula[1],this->field->BSDF);
  if (parameters.field_formula[2])
    this->field->add_field_formulaSD('z', parameters.field_formula[2],this->field->BSDF);

  }
  else{
  if (parameters.field_formula[0])
    this->field->add_field_formula('x', parameters.field_formula[0]);
  if (parameters.field_formula[1])
    this->field->add_field_formula('y', parameters.field_formula[1]);
  if (parameters.field_formula[2])
    this->field->add_field_formula('z', parameters.field_formula[2]);
  }
  

  //std::cout<<"past Initial pulse"<<std::endl;

  // Diffusion
  this->boundary->SetDiffusion(parameters.diffusion);
  // Gravity
  if (parameters.gravity != 0.) {
    this->boundary->SetGravity(parameters.gravity);
  }

  // Initialize factory
  
  //delete factory; deleting factory doesn't work. but it is a memory leak, I think.
  this->factory = new ParticleFactory();
  
  this->factory->parameters = parameters;
  this->factory->field = field;
  this->factory->boundary = boundary;
  this->factory->scattering = scattering;
  
  for (int i=0; i < 4; i++) {
    factory->vDistribution[i] = vDistribution[i];
  }
  //std::cout<<"Done reloading paraemters"<<std::endl;
}

void Run::useQuiteBoltzmannDistribution(Double_t temperature) {

  TF1 *f1;

  const Double_t boltzmann = 1.38066e-23; // J K^-1
  const Double_t MHe3 = 2.2 * 3.016029 * 1.66054e-27; // kg !!effective Mass (M*=2.2M)

  //cout << "Using Boltzmann (T = " << temperature << " K)" << endl;

  Double_t temp = MHe3 / (2. * temperature * boltzmann);

  f1 = new TF1("boltzmann", "(x^2) * exp(-(x^2) * [0])", 0, 
         10. * sqrt(1. / temp));
  // limit set to 10 v_p (the most probable velocity)
  // 4.2e-43 integrated probability beyond that velocity

  f1->SetNpx(2000);
  
  f1->SetParameter(0,temp);

  this->vDistribution[0] = f1;


  temp = sqrt(8. * boltzmann * temperature / (MHe3 * 3.14159));

  this->mfp = 3. * (1.6e-4 / pow(temperature, 7)) / temp;

  temp = sqrt(3. * boltzmann * temperature / MHe3);

  this->tau_C = 3. * (1.6e-4 / pow(temperature, 7)) / (temp*temp);

//  cout << "MFP = " << this->mfp << " m" << endl;

  this->scattering = new Scattering(1);
  this->scattering->Set(this->tau_C);

 // cout << "Scattering class enabled with \\tau_C = " << this->tau_C << endl;

  return;


}


