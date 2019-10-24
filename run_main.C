#include "Run.h"

#include <iostream>
#include <stdlib.h>

extern int RUN_EXIT_ASAP;
extern int RUN_NEEDS_SAVE;

int main (int argc, char **argv) {

  Run *r;

  // Syntax:
  // run <runID> [OPTIONS]

  int options[] = {0,0,-1,0,-1};
  // 0 -p
  // 1 -c
  // 2 -n
  // 3 -novc
  // 4 -numerical

  int runfileExists = 0;

  char *par_file;

  // Get run ID from argument line
  if (argc < 2) {
    cout << "Please enter run ID" << endl;
    return 0;
  }
  r = new Run(atoi(argv[1]));

  Run::parse_argument_line(argc, argv, options, &par_file);

  // Check if the run already exists
  runfileExists = Run::FileExists(r->filename);

  // New simulation starts
  if (runfileExists == 0) {
    // option -p needed. -c not wanted
    if (options[0] == 0 || options[1] != 0) {
      cout << "Simulation does not exist please use the syntax:"
	   << endl << argv[0] << "runID -p param_file" << endl
	   << "\t(optional arg:\t-n nEntries)" << endl
	   << "\t(optional arg:\t-numerical method_id)" << endl
	   << "Do not include other arguments." << endl;
      return 1;
    }

    // get run parameters
    r->parameters.ParseParameterFile(par_file);
    if (options[2] > 0) {
      cout << "-n option overriding nEntries. From "
	   << r->parameters.nEntries << " to " << options[2] << "." << endl;
      r->parameters.nEntries = options[2];
    }
    if (r->parameters.VerifyParameters() > 0) {
      cout << "Please revise parameter file" << endl;
      return 1;
    }
    cout << "// Parameters\n///////////////////////" << endl;
    r->parameters.PrintParameters();
    cout << "/////////////\n// End Parameters" << endl;
    r->nCurrent = 0;
  }
  else { // need to modify to allow continuation of runs
    if (options[1] != 1) {
      cout << "Need -c parameter to continue run" << endl;
      return 0;
    }
    cout << "Run already exists" << endl;
    if (r->loadRun() != 0) {
      cout << "Could not load Run" << endl;
      return -1;
    }
    cout << r->nCurrent << " entries" << endl;
    if (options[2] > 0) {
      cout << "Can continue to " << options[2] << endl;
      r->parameters.nEntries = options[2];
    }
    else {
      cout << "Use -c argument to continue the run." << endl;
      return 1;
    }

    r->parameters.PrintParameters();
    //return 1;
  }
  // Numerical method
  if (options[4] >= 0) {
    r->parameters.numerical_method = options[4];
  }

  // Initialize run
  r->InitializeRun();
  r->InitializeTFile();
  set_signal_handles();

  int not_written = 0;
  cout << "Starting run..." << endl;
  extern int RUN_EXIT_ASAP;
  extern int RUN_NEEDS_SAVE;
  for (int i=r->nCurrent; i<r->parameters.nEntries && !RUN_EXIT_ASAP; i++) {

    if(r->runNeutron() == 0) {
      cout << "." << flush;
      not_written++;
      RUN_NEEDS_SAVE = 1;
      if (not_written > 1000) {
	r->file->Write(NULL, TObject::kWriteDelete);
	RUN_NEEDS_SAVE = 0;
	not_written = 0;
      }
    }
    else {
      cout << "discarded one neutron" << endl;
      // If this ever happens, it would be necessary to
      // restore the random number generator
      i--;
    }

  }

  if (not_written > 0) {
    r->file->Write(NULL, TObject::kWriteDelete);
  }

  delete r;

  cout << endl;

  cout << "Exit ok" << endl;
}
