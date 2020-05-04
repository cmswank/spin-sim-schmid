#!/usr/bin/env python

import os
import sys
import getopt
import numpy as np

from ROOT import gSystem
gSystem.Load(os.path.dirname(__file__)+'/runShared')

from ROOT import Run, RunParameters

def create_run(run_id, opts):
    run = Run(run_id)
    if Run.FileExists(run.filename):
        print "%s already exists" % (run.filename, )
        if '-c' not in opts:
            raise ValueError('Need -c parameter to continue existing run')
        if '-p' in opts:
            print "W: ignoring parameter file %s" % (opts['-p'],)
        if run.loadRun() != 0:
            raise ValueError("Could not load run %d" % (run_id,))
        print "%d entries found" % (run.nCurrent,)
        if '-n' in opts:
            n = int(opts['-n'])
            if n <= run.nCurrent:
                raise ValueError('use -n to continue to a number larger than %d' % (run.nCurrent,))
            run.parameters.nEntries = n
        print "continuing to %d" % (run.parameters.nEntries,)
    else:
        if '-c' in opts:
            raise ValueError('-c fails. Nothing to continue, or %s not found' % (run.filename,))
        if '-p' not in opts:
            raise ValueError('need to specify parameter file with -p')
        run.parameters.ParseParameterFile(opts['-p'])
        if '-n' in opts:
            n = int(opts['-n'])
            print '-n used. Overriding nEntries %d -> %d' % (run.parameters.nEntries, n)
            run.parameters.nEntries = n
        run.nCurrent = 0
    run.parameters.PrintParameters()
    if '--numerical' in opts:
        run.parameters.numerical_method = int(opts['--numerical'])
    run.InitializeRun()
    run.InitializeTFile()
    return run

def main_loop(run):
    print "Starting run %d..." % (run.runID,)
    while run.nCurrent < run.parameters.nEntries:
        rn=run.nCurrent
        if rn >=0: #careful. 

            #run.parameters.Pulseparam[1]=rn//2+1;
            ##E field Variation
            if rn%2==0:
                run.parameters.E0[0]=-run.parameters.E0[0];
                if rn ==4:
                    run.parameters.E0[0]=0;

            run.parameters.simType=rn%2
            if rn%2==0:
                run.parameters.spinSettings[0]=0.
                run.parameters.spinSettings[1]=np.cos(0.48/2)
                run.parameters.spinSettings[2]=np.sin(0.48/2)
                run.parameters.edm=1E-14;
                run.parameters.speed=5;
            else:
                run.parameters.spinSettings[0]=0.
                run.parameters.spinSettings[1]=np.cos(0.48/2)
                run.parameters.spinSettings[2]=-np.sin(0.48/2)
                run.parameters.edm=0;
            run.parameters.seed=678742112
            run.ReloadParameters()
        ret = run.runNeutron()

        if ret == 0:
            sys.stdout.write(str(rn)+' ')
            sys.stdout.flush()
        else:
            print "discarded one neutron"

if __name__ == '__main__':
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], 'cp:n:', ['novc', 'numerical='])
    except getopt.GetoptError as e:
        print e
        sys.exit(1)
    opts = dict(opts)
    try:
        run_id = int(args[0])
    except IndexError:
        print "Please enter run ID"
        sys.exit(1)
    except ValueError:
        print "Please enter valid integer run ID"
        sys.exit(1)

    try:
        run = create_run(run_id, opts)
    except Exception as e:
        print "There was an error creating the run"
        print e
        sys.exit(1)

    try:
        main_loop(run)
    except KeyboardInterrupt:
        print "^C caught. Letting file finish saving..."
    else:
        print "\ndone"
    finally:
        run.WriteTFile()
        print "exit ok"
        del run

    sys.exit(0)
