#!/usr/bin/env python

import os
import sys
import getopt
import numpy as np

from ROOT import gSystem
gSystem.Load('/data1/cmswank/spin-sim-xliu/runShared.so')

from ROOT import Run, RunParameters

def create_run(run_id, opts):
    run = Run(run_id)
    if Run.FileExists(run.filename):
        print("%s already exists" % (run.filename, ))
        if '-c' not in opts:
            raise ValueError('Need -c parameter to continue existing run')
        if '-p' in opts:
            print("W: ignoring parameter file %s" % (opts['-p'],))
        if run.loadRun() != 0:
            raise ValueError("Could not load run %d" % (run_id,))
        print("%d entries found" % (run.nCurrent,))
        if '-n' in opts:
            n = int(opts['-n'])
            if n <= run.nCurrent:
                raise ValueError('use -n to continue to a number larger than %d' % (run.nCurrent,))
            run.parameters.nEntries = n
        print("continuing to %d" % (run.parameters.nEntries,))
    else:
        if '-c' in opts:
            raise ValueError('-c fails. Nothing to continue, or %s not found' % (run.filename,))
        if '-p' not in opts:
            raise ValueError('need to specify parameter file with -p')
        run.parameters.ParseParameterFile(opts['-p'])
        if '-n' in opts:
            n = int(opts['-n'])
            print('-n used. Overriding nEntries %d -> %d' % (run.parameters.nEntries, n))
            run.parameters.nEntries = n
        run.nCurrent = 0
    run.parameters.PrintParameters()
    if '--numerical' in opts:
        run.parameters.numerical_method = int(opts['--numerical'])
    run.InitializeRun()
    run.InitializeTFile()
    return run

def main_loop(run):
    print("Starting run %d..." % (run.runID,))
    while run.nCurrent < run.parameters.nEntries:
        rn=run.nCurrent
        rnf=run.nCurrent*1.0000
        norsf=run.parameters.nEntries*1.0000;
        if rn >=0: #careful. 

            
            #B1=0.000250021223525*(1+5E-5*rn//2)
                #print "Bz1 = "+str(B1)+","
            #run.EraseFieldFormula();
            #run.parameters.field_assignment("BzAdd",str(B1))
            #run.parameters.Pulseparam[1]=rn//2+1;
            ##E field Variation
            #if rn%2==0:
                #run.parameters.E0[0]=-run.parameters.E0[0];
                #run.parameters.E0[0]=0.0;
                

            #run.parameters.E0[0]=0.;        
            run.parameters.simType=rn%2
            if rn%2==0:
                #run.parameters.seed=282*(rn+1)    
                #4.02497057364e-05; best critically dressed for 5.2 uT. 103 Hz. 8 significant digits.
                #Best for robust dressing at 5.2 uT 200 Hz mod, B1=4.205e-5; modulation amplitude = 3.497!. 
                #Best for robust dressing at 5.2 uT 300 Hz mod, B1=4.42958e-5; modulation amplitude = 3.32973!. 
                #best for robust dresssing at 5.2 uT 550Hz mod, B1=5.3491388e-05; mod amp= 2.591036

                #B1=4.205e-5#5.3491388e-05#4.205e-5#5.3491388E-5#4.205e-5#5.347e-05*(1+1E-3*(rn/2-run.parameters.nEntries/4)/(run.parameters.nEntries/2));#4.02444353976e-05*(1+1E-4*(rn/2-run.parameters.nEntries/4)/(run.parameters.nEntries/2));#0.000250021223525*(1+1E-3*(rn/2-run.parameters.nEntries/4)/(run.parameters.nEntries/2))#
                #print "Bz1 = "+str(B1)+"\n"
                #run.EraseFieldFormula();
                #run.parameters.field_assignment("BzAdd",str(B1))
                
                #He3 
                #run.parameters.SDparam[4]=3.57E-5/B1*3.497 #*(1+1E-3*(rn/2-run.parameters.nEntries/4)/(run.parameters.nEntries/2));
                #print "x1 = "+str(run.parameters.SDparam[4])+"\n"
                #print "Mod amp = "+str(run.parameters.SDparam[4]*B1/3.57E-5)+"\n"
                #6 7 Bz1 = 4.1915874e-05

                #
                #run.parameters.Pulseparam[3]=1.245#*(1+.05*(rn/2-run.parameters.nEntries/4)/(run.parameters.nEntries/2))
                #run.parameters.Pulseparam[2]=0.708#*(1+0.1*(rn/2-run.parameters.nEntries/4)/(run.parameters.nEntries/2))#7 is noise rng seed. 
                #Neutrons:
                #3845.309408, 6283.18530718
                #
                #phase=np.pi/2;#np.pi/4;#np.pi/2;#np.pi/4.0;
                #phi=0;

                phase=np.pi/8;#np.pi/2;#np.pi/4.0;
                phi=0.;
                #run.parameters.B0[0]=5.15E-6*(1+1E-2*(x[1]-1));
                #run.parameters.spinSettings[0]=0.
                #run.parameters.spinSettings[1]=np.sqrt(1.-np.square(run.parameters.spinSettings[0]))*np.sin(phi+phase)#+phase/2.0)
                #run.parameters.spinSettings[2]=np.sqrt(1.-np.square(run.parameters.spinSettings[0]))*np.cos(phi+phase)#+phase/2.0)
                run.parameters.spinSettings[0]=-0.001441682659922;
                run.parameters.spinSettings[1]=-0.405022553486012;  
                run.parameters.spinSettings[2]=0.914305557633103;

                #run.parameters.simEDM=0#1E-23 #this is e meters. (ugghhh)
                #run.parameters.edm=0#1E-23;
                run.parameters.speed=5
                #if rn%4==0:
                    #run.parameters.SDparam[2]=np.pi*(rnf/norsf)
                    #run.parameters.SDparam[2]=0.0;
                
                #run.parameters.Pulseparam[7]=82*(rn+1) #7 is noise rng seed.
                #print "B1phi = " +str(run.parameters.SDparam[2])+"\n"
                #run.parameters.Pulseparam[1]=(rn//2+2)  
                #offset in another cell, 
                
            else:
                #3He:
                #run.parameters.spinSettings[0]=0;#np.cos(np.pi/4);#0.
                #run.parameters.spinSettings[1]=1;#np.sin(np.pi/4);#sin(phi-phase/2.0)#np.pi/4)
                #run.parameters.spinSettings[2]=0.;#np.cos(phi+phase/2.0)#np.pi/4)
                phi=0.;
                phase=np.pi/8
                #run.parameters.B0[0]=5.15E-6*(1+1E-2*(x[1]-1));
                #run.parameters.spinSettings[0]=0.
                #run.parameters.spinSettings[1]=np.sqrt(1.-np.square(run.parameters.spinSettings[0]))*np.sin(phi-phase)#+phase/2.0)
                #run.parameters.spinSettings[2]=np.sqrt(1.-np.square(run.parameters.spinSettings[0]))*np.cos(phi-phase)#+phase/2.0)
                run.parameters.spinSettings[0]=0.001616473758831;
                run.parameters.spinSettings[1]=0.360218604798280;
                run.parameters.spinSettings[2]=0.932866519803218;


                run.parameters.edm=0.0
                run.parameters.simEDM=0.0
                #run.parameters.E
                #print "phi = " +str(run.parameters.SDparam[2])+"\n"
                #7 is noise rng seed. 
                #run.parameters.Pulseparam[7]=82*(rn+1) #7 is noise rng seed. 
             #678742112  #seed is trajectory rng. 
            run.ReloadParameters();
            
        ret = run.runNeutron()


        if ret == 0:
            sys.stdout.write(str(rn)+' ')
            sys.stdout.flush()
        else:
            print("discarded one neutron")

if __name__ == '__main__':
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], 'cp:n:', ['novc', 'numerical='])
    except getopt.GetoptError as e:
        print(e)
        sys.exit(1)
    opts = dict(opts)
    try:
        run_id = int(args[0])
    except IndexError:
        print("Please enter run ID")
        sys.exit(1)
    except ValueError:
        print("Please enter valid integer run ID")
        sys.exit(1)

    try:
        run = create_run(run_id, opts)
    except Exception as e:
        print("There was an error creating the run")
        print(e)
        sys.exit(1)

    try:
        main_loop(run)
    except KeyboardInterrupt:
        print("^C caught. Letting file finish saving...")
    else:
        print("\ndone")
    finally:
        run.WriteTFile()
        print("exit ok")
        del run

    sys.exit(0)
