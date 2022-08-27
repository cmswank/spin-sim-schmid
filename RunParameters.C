#include "RunParameters.h"

#include "Run.h"
#include "Vector.h"

RunParameters::RunParameters() : tree(0), simType(-1), simGamma(0),
				 simEDM(0), seed(0), numerical_method(5),
				 max_error(1e-9), max_angle(PI/8),
				 totalTime(0), offsetTime(0),
				 nBins(0), randBins(0),
				 speed(0), temperature(0),
				 geometry(-1), radius(-1), diffusion(0),
				 gravity(0), uniform_gradient(0),
				 simple_gradient(0), simple_quad_gradient(0),
				 vdistrexpo(-1), no_normalize(0) {
  version[0] = RUN_VERSION[0];
  version[1] = RUN_VERSION[1];
  version[2] = RUN_VERSION[2];
  boxLow[0] = boxLow[1] = boxLow[2] = 0;
  boxHigh[0] = boxHigh[1] = boxHigh[2] = 0;
  spinSettings[0] = spinSettings[1] = spinSettings[2] = -10;
	
  SDparam[0] = SDparam[1] = 0;
  B0[0] = B0[1] = B0[2] = 0;
  E0[0] = E0[1] = E0[2] = 0;

  // other parameters
  field_formula[0] = field_formula[1] = field_formula[2] = 0;
}

RunParameters::~RunParameters() {
  if (uniform_gradient)
    delete[] uniform_gradient;
  if (simple_gradient)
    delete[] simple_gradient;
  if (simple_quad_gradient)
    delete[] simple_quad_gradient;
}

void RunParameters::BranchTTree(TTree *tree) {
  this->tree = tree;
  tree->Branch("par1", store_par1, "x/D:y:z");
  tree->Branch("id", &store_par_id, "id/I");
  tree->Branch("info", store_par_info, "info/C");
}

void RunParameters::SetTTreeBranches(TTree *tree) {
  this->tree = tree;
  tree->SetBranchAddress("par1", store_par1);
  tree->SetBranchAddress("id", &store_par_id);
  tree->SetBranchAddress("info", store_par_info);
}

void RunParameters::PrintParameters() {
  cout << "// ----- Simulation" << endl;
  cout << "simType = " << simType << ";" << endl;
  cout << "simGamma = " << simGamma << ";" << endl;
  cout << "seed = " << seed << ";" << endl;
  cout << "runTime = " << totalTime << ";" << endl;
  if (offsetTime > 0.) {
    cout << "offsetTime = " << offsetTime << ";" << endl;
  }
  cout << "nEntries = " << nEntries << ";" << endl;
  cout << "nBins = " << nBins << ";" << endl;
  cout << "randBins = " << randBins << ";" << endl;
  cout << "numerical = " << numerical_method << ";" << endl;
  cout << "maxError = " << max_error << ";" << endl;
  cout << endl;

  cout << "// ----- Velocity" << endl;
  cout << "velocity = " << speed << ";" << endl;
  cout << "temperature = " << temperature << ";" << endl;
  if (vdistrexpo >= 0.) {
    cout << "vdistrexpo = " << vdistrexpo << ";" << endl;
  }
  cout << endl;

  cout << "// ----- spinSettings" << endl;
  if (spinSettings[1] < -1.1) {
    cout << "spinSet = " << spinSettings[0] << ";" << endl;
  }
  else {
    cout << "spinSet = " << spinSettings[0] << ", " << spinSettings[1] << ", " << spinSettings[2] << ";" << endl;
  }
  cout << endl;

  cout << "// ----- Geometry" << endl;
  cout << "geometry = " << geometry << ";" << endl;

  if (geometry == 1 || geometry == 3 || geometry == 4) {
    cout << "boxLow = " << boxLow[0] << ", " << boxLow[1] << ", " << boxLow[2] << ";" << endl;
    cout << "boxHigh = " << boxHigh[0] << ", " << boxHigh[1] << ", " << boxHigh[2] << ";" << endl;
  }
  if (geometry == 0 || geometry == 4) {
    cout << "radius = " << radius << ";" << endl;
  }

  cout << "diffusion = " << diffusion << ";" << endl;
  if (gravity != 0.) {
    cout << "gravity = " << gravity << ";" << endl;
  }
  cout << endl;

  cout << "// ----- Fields" << endl;

  cout << "B0 = " << B0[0] << ", " << B0[1] << ", " << B0[2] << ";" << endl;

  cout << "E0 = " << E0[0] << ", " << E0[1] << ", " << E0[2] << ";" << endl;
  if (uniform_gradient) {
    cout << "UniformG = "
	 << uniform_gradient[0] << ", "
	 << uniform_gradient[1] << ", "
	 << uniform_gradient[2] << ", "
	 << uniform_gradient[3] << ", "
	 << uniform_gradient[4] << ", "
	 << uniform_gradient[5] << ", "
	 << uniform_gradient[6] << ", "
	 << uniform_gradient[7] << ", "
	 << uniform_gradient[8] << ";" << endl;
  }
  cout << endl;

  if (field_formula[0])
    cout << "Bx += " << field_formula[0]->GetExpFormula() << endl;
  if (field_formula[1])
    cout << "By += " << field_formula[1]->GetExpFormula() << endl;
  if (field_formula[2])
    cout << "Bz += " << field_formula[2]->GetExpFormula() << endl;

  cout << endl;

  cout << "// -----\n";
}

void RunParameters::LoadParameters() {
  if (!tree) {
    cout << "W: No tree to load parameters from" << endl;
    return;
  }
  SetTTreeBranches(tree);

  tree->GetEntry(0);
  version[0] = (int)(store_par1[0]);
  version[1] = (int)(store_par1[1]);
  version[2] = (int)(store_par1[2]);

  tree->GetEntry(1);
  simType = (Int_t)(store_par1[0]);
  seed = (Int_t)(store_par1[1]);
  simGamma = store_par1[2];

  tree->GetEntry(2);
  nBins = (Int_t)(store_par1[0]);
  totalTime = store_par1[1];
  randBins = (Int_t)(store_par1[2]);

  tree->GetEntry(3);
  B0[0] = store_par1[0];
  B0[1] = store_par1[1];
  B0[2] = store_par1[2];

  tree->GetEntry(4);
  E0[0] = store_par1[0];
  E0[1] = store_par1[1];
  E0[2] = store_par1[2];

  tree->GetEntry(5);
  geometry = (Int_t)(store_par1[0]);
  radius = store_par1[1];
  simEDM = store_par1[2];

  tree->GetEntry(6);
  boxLow[0] = store_par1[0];
  boxLow[1] = store_par1[1];
  boxLow[2] = store_par1[2];

  tree->GetEntry(7);
  boxHigh[0] = store_par1[0];
  boxHigh[1] = store_par1[1];
  boxHigh[2] = store_par1[2];

  tree->GetEntry(8);
  diffusion = store_par1[0];
  gravity = store_par1[1];

  tree->GetEntry(9);
  speed = store_par1[0];
  temperature = store_par1[1];
  offsetTime = store_par1[1];

  tree->GetEntry(10);
  spinSettings[0] = store_par1[0];
  spinSettings[1] = store_par1[1];
  spinSettings[2] = store_par1[2];

  tree->GetEntry(11);
  numerical_method = (Int_t)(store_par1[0]);
  max_error = store_par1[1];
  max_angle = store_par1[2];

  uniform_gradient = new Double_t[9]();
  tree->GetEntry(12);
  uniform_gradient[0] = store_par1[0];
  uniform_gradient[1] = store_par1[1];
  uniform_gradient[2] = store_par1[2];

  tree->GetEntry(13);
  uniform_gradient[3] = store_par1[0];
  uniform_gradient[4] = store_par1[1];
  uniform_gradient[5] = store_par1[2];

  tree->GetEntry(14);
  uniform_gradient[6] = store_par1[0];
  uniform_gradient[7] = store_par1[1];
  uniform_gradient[8] = store_par1[2];
  // if uniform gradient is null leave it 0
  int i;
  for (i=0; i<9; i++) {
    if (uniform_gradient[i] != 0)
      break;
  }
  if (i == 9) {
    delete[] uniform_gradient;
    uniform_gradient = 0;
  }
}

void RunParameters::SaveParameters() {
  // create ttree if it does not exist
  if (!tree) {
    tree = new TTree("par", "Parameter container");
    BranchTTree(tree);
  }

  SaveParVector(0, (Double_t)version[0], (Double_t)version[1], (Double_t)version[2], "Version");

  SaveParVector(1, (Double_t)simType, (Double_t)seed, simGamma, "simType, seed, gamma");

  SaveParVector(2, (Double_t)nBins, totalTime, (Double_t)randBins, "nBins, t, randBins");

  SaveParVector(3, B0, "Magnetic Field (T)");
  SaveParVector(4, E0, "Electric Field (V/m)");

  SaveParVector(5, (Double_t)geometry, radius, simEDM, "Geometry, radius, simEDM");

  SaveParVector(6, boxLow, "Box Low Limits");
  SaveParVector(7, boxHigh, "Box High Limits");

  SaveParVector(8, diffusion, gravity, 0, "diffusion, gravity");

  SaveParVector(9, speed, temperature, offsetTime, "Speed, Temperature");

  SaveParVector(10, spinSettings, "Spin Settings");

  SaveParVector(11, (Double_t)numerical_method, max_error, max_angle, "Numerical method, maxError, maxAngle");

  if (uniform_gradient) {
    SaveParVector(12, uniform_gradient, "dBx");
    SaveParVector(13, uniform_gradient+3, "dBy");
    SaveParVector(14, uniform_gradient+6, "dBz");
  }
  else {
    SaveParVector(12, 0,0,0, "dBx");
    SaveParVector(13, 0,0,0, "dBy");
    SaveParVector(14, 0,0,0, "dBz");
  }
}

void RunParameters::SaveParVector(Int_t id, Double_t x, Double_t y, Double_t z, const Char_t *description) {
  store_par1[0] = x;
  store_par1[1] = y;
  store_par1[2] = z;

  if (description)
    sprintf(store_par_info, "%s", description);
  else
    sprintf(store_par_info, "N/A");

  tree->Fill();
}

void RunParameters::SaveParVector(Int_t id, Double_t *vals, const Char_t *description) {
  SaveParVector(id, vals[0], vals[1], vals[2], description);
}

int RunParameters::check_assignment_length(Reading *param, Int_t args) {
  if (param->getVarLength() < args) {
    cout << "ERROR: " << param->getVarName() << " needs "
	 << args << " args\n";
    return -1;
  }
  if (param->getVarLength() > args) {
    cout << "Warning: " << param->getVarName() << " only needs "
	 << args << " args\n";
    return 1;
  }
  return 0;
}

/**
 * Assign parameter values to class variables
 * @param Reading *param : the Reading class parameter in question
 * @param Int_t args : the required number of arguments
 * @param Double_t *var : the class variable (or array)
 */
int RunParameters::assignment(Reading *param, Int_t args, Double_t *var) {
  if (check_assignment_length(param, args) < 0)
    return 0;
  for (int i=0; i<args; i++) {
    var[i] =  param->getVarVector()[i];
  }
  return 1;
}

int RunParameters::assignment(Reading *param, Int_t args, Int_t *var) {
  if (check_assignment_length(param, args) < 0)
    return 0;
  for (int i=0; i<args; i++) {
    var[i] =  (Int_t)(param->getVarVector()[i] + 0.1);
  }
  return 1;
}

int RunParameters::assignment(Reading *param, Int_t args, UInt_t *var) {
  if (check_assignment_length(param, args) < 0)
    return 0;
  for (int i=0; i<args; i++) {
    var[i] =  (UInt_t)(param->getVarVector()[i] + 0.1);
  }
  return 1;
}

int RunParameters::assignment(std::string *expression, Int_t args, Char_t *var) {
  //if (check_assignment_length(expression, args) < 0)
   // return 0;
	
 var=(Char_t*)expression->c_str();
 std::cout<<var[0]<<"\n";
  //for (int i=0; i<args; i++) {
//	  std::cout<<param->getVarVector()[i]<<"\n";
 //   var[i] =  (Char_t)(param->getVarVector()[i]);
  //}
  return 1;
}


Int_t RunParameters::field_assignment (const char *field_component, const char *expression) {

  TFormula *field_function;

  field_function = new TFormula("temp",expression);
  if (field_function->Compile() == 1) {
    cout << "Error in function " << field_component
	 << " = \"" << expression << "\"." << endl
	 << "Function will not be used." << endl;
    delete field_function;
    return 0;
  }

  TFormula **p_field_temp = 0;
  TString temp_string, function_name;
  if (strcmp(field_component,"BxAdd") == 0) {
    p_field_temp = &this->field_formula[0];
    function_name = "field_formula_x";
  }
  else if (strcmp(field_component,"ByAdd") == 0) {
    p_field_temp = &this->field_formula[1];
    function_name = "field_formula_y";
  }
  else if (strcmp(field_component,"BzAdd") == 0) {
    p_field_temp = &this->field_formula[2];
    function_name = "field_formula_z";
  }


  if (!(*p_field_temp)) {
    (*p_field_temp) = new TFormula(function_name,
				   field_function->GetExpFormula());
  }
  else {
    temp_string = (*p_field_temp)->GetExpFormula();
    temp_string.Append(" + ");
    temp_string.Append(field_function->GetExpFormula());
    delete (*p_field_temp);
    (*p_field_temp) = new TFormula(function_name, temp_string);
  }


  return 1;
}

int RunParameters::ParseParameterFile (const char *filename) {
  Reading *param = new Reading();

  char *tempStr, *expression;
  int varType;

  if (param->open(filename) != 0) {
    cout << "Error opening parameter file\n";
    return -1;
  }

  cout << "\nStarting to parse file " << filename << endl;

  while ((varType = param->getNextAssignment())) {
    tempStr = param->getVarName();
    if (varType == 1) {

      // Geometry
      if (strcmp(tempStr,"geometry") == 0)
	assignment(param, 1, &geometry);

      else if (strcmp(tempStr,"radius") == 0)
	assignment(param, 1, &radius);

      else if (strcmp(tempStr,"boxLow") == 0)
	assignment(param, 3, boxLow);

      else if (strcmp(tempStr,"boxHigh") == 0)
	assignment(param,3,boxHigh);

      else if (strcmp(tempStr,"diffusion") == 0)
	assignment(param,1,&diffusion);

      else if (strcmp(tempStr,"gravity") == 0)
	assignment(param,1,&gravity);

      // Fields
      else if (strcmp(tempStr,"B0") == 0)
	assignment(param,3,B0);

      else if (strcmp(tempStr,"UniformG") == 0) {
	uniform_gradient = new Double_t[9]();
	assignment(param,9,uniform_gradient);
      }

      else if (strcmp(tempStr,"SimpleGrad") == 0) {
	simple_gradient = new Double_t[3];
	assignment(param,3,simple_gradient);
      }

      else if (strcmp(tempStr,"SimpleQuadGrad") == 0) {
	simple_quad_gradient = new Double_t[3];
	assignment(param,3,simple_quad_gradient);
      }

      else if (strcmp(tempStr,"E0") == 0)
	assignment(param,3,E0);

      // Simulation
      else if (strcmp(tempStr, "simType") == 0)
	assignment(param,1,&simType);

      else if (strcmp(tempStr, "simGamma") == 0)
	assignment(param,1,&simGamma);

      else if (strcmp(tempStr, "edm") == 0)
	assignment(param,1,&simEDM);

      else if (strcmp(tempStr, "noNormalize") == 0)
	assignment(param,1,&no_normalize);

      else if (strcmp(tempStr,"nEntries") == 0)
	assignment(param,1,&nEntries);

      else if (strcmp(tempStr,"nBins") == 0)
	assignment(param,1,&nBins);

      else if (strcmp(tempStr,"randBins") == 0)
	assignment(param,1,&randBins);

      else if (strcmp(tempStr, "seed") == 0)
	assignment(param,1,&seed);

      else if (strcmp(tempStr, "runTime") == 0)
	assignment(param,1,&totalTime);

      else if (strcmp(tempStr, "offsetTime") == 0)
	assignment(param,1,&offsetTime);

      else if (strcmp(tempStr, "numerical") == 0)
	assignment(param,1,&numerical_method);

      else if (strcmp(tempStr, "maxError") == 0)
	assignment(param,1,&max_error);

      else if (strcmp(tempStr, "maxAngle") == 0)
	assignment(param,1,&max_angle);

      // Velocity
      else if (strcmp(tempStr, "velocity") == 0)
	assignment(param,1,&speed);

      else if (strcmp(tempStr, "temperature") == 0)
	assignment(param,1,&temperature);

      else if (strcmp(tempStr, "vdistrexpo") == 0)
	assignment(param,1,&vdistrexpo);
	
	//Spin Dressing added by C.Swank	
	else if (strcmp(tempStr,"SpinDressing")==0)
	  assignment(param,8,SDparam);

    //Interp added by C.Swank  
  else if (strcmp(tempStr,"FieldInterp")==0)
    assignment(param,10,Interpparam);

	
	//Initial pulse added by E. Webb
  else if (strcmp(tempStr,"InitialPulse")==0)
    assignment(param,14,Pulseparam);

  else if (strcmp(tempStr,"noise")==0)
    assignment(param,1,&noise);
      
		// Spin
      else if (strcmp(tempStr, "spinSet") == 0) {
	if (param->getVarLength() == 3) {
	  assignment(param,3,spinSettings);
	}
	else {
	  assignment(param,1,spinSettings);
	}
      }


      // not recognized
      else {
	cout << "Error: parameter <" << tempStr << "> not recognized\n";
	cout << "---- " << param->getLineInput() << endl;
      }

    } // end if (varType == 1)
    else if (varType == 2) {
      expression = param->getExpression();

      if (strcmp(tempStr, "BxAdd") == 0) {
	field_assignment("BxAdd", expression);
      }
      else if (strcmp(tempStr, "ByAdd") == 0) {
	field_assignment("ByAdd", expression);
      }
      else if (strcmp(tempStr, "BzAdd") == 0) {
	field_assignment("BzAdd", expression);
      }
		
	  else if (strcmp(tempStr,"FieldIntDir")==0){
	  InterpDir = (std::string)param->getExpression();
	  //assignment(expres,3,InterpDir);
	  }
   	
      // not recognized
      else {
	cout << "Error: parameter <" << tempStr << "> not recognized\n";
	cout << "---- " << param->getLineInput() << endl;
      }
    }
  }

  // deconstruct the reader
  delete param;

  return 0;
}

// Verify that the parameters are good for the simulation
int RunParameters::VerifyParameters() {
  int returnValue = 0;

  if (simType == 0) {
    if (speed < 0) {
      cout << "ERROR: speed not set" << endl;
      returnValue++;
    }
  }
  else if (simType == 1) {
    if (temperature < 0) {
      cout << "ERROR: temperature not set" << endl;
      returnValue++;
    }
  }

  // Check numerical method
  switch (numerical_method) {
  case 0:
  case 1:
  case 5:
    break;
  default:
    cout << "ERROR: Numerical method " << numerical_method
	 << " is not defined." << endl;
    returnValue++;
  }

  // Check boundary
  if ((geometry == 0 || geometry == 2 || geometry == 4) && radius < 0) {
    cout << "ERROR: radius not set" << endl;
    returnValue++;
  }
  else if ((geometry == 1) &&
	   (boxLow[0] > boxHigh[0] ||
	    boxLow[1] > boxHigh[1] ||
	    boxLow[2] > boxHigh[2])) {
    cout << "ERROR: box not properly set" << endl;
    returnValue++;
  }
  else if ((geometry == 3) &&
	   (boxLow[0] > boxHigh[0] ||
	    boxLow[1] > boxHigh[1])) {
    cout << "ERROR: 2D-box not properly set" << endl;
    returnValue++;
  }

  // Check Fields
  if (B0[0] == 0. && B0[1] == 0. && B0[2] == 0.) {
    cout << "Warning: B field is zero" << endl;
  }

  if (E0[0] == 0. && E0[1] == 0. && E0[2] == 0.) {
    cout << "Warning: E field is zero" << endl;
  }

  // Simulation time
  if (totalTime < 0) {
    cout <<"ERROR: simulation time not set" << endl;
    returnValue++;
  }

  // nEntries
  if (nEntries < 0) {
    cout << "ERROR: nEntries not defined, can use -n argument." << endl;
    returnValue++;
  }

  // Spin Settings
  if (spinSettings[0] < -1.1) {
    cout << "ERROR: No spin set." << endl;
    returnValue++;
  }

  return returnValue;

}
