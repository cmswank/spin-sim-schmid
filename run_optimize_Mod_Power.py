#!/usr/bin/env python

import os
import sys
import getopt
import numpy as np
from scipy.optimize import minimize

from ROOT import gSystem

gSystem.Load('/data1/cmswank/spin-sim-xliu/runShared.so') #the fancy stuff broke after c libraries updated. 

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

################################## Optimization Loop  !!!!!!!!!!

def single_loop(x1,run):
    
    print("Starting run %d for x1 = %f..." % (run.runID,x1[0]))
    Sxn=np.empty([run.parameters.nEntries//2,run.parameters.nBins+1]);
    Syn=np.empty([run.parameters.nEntries//2,run.parameters.nBins+1]);
    Szn=np.empty([run.parameters.nEntries//2,run.parameters.nBins+1]);
    Sx3=np.empty([run.parameters.nEntries//2,run.parameters.nBins+1]);
    Sy3=np.empty([run.parameters.nEntries//2,run.parameters.nBins+1]);
    Sz3=np.empty([run.parameters.nEntries//2,run.parameters.nBins+1]);
    tsim=np.empty([run.parameters.nEntries//2,run.parameters.nBins+1]);
    run.nCurrent=0;


    while run.nCurrent < run.parameters.nEntries:
        rn=run.nCurrent
        rnf=run.nCurrent*1.0000
        norsf=run.parameters.nEntries*1.0000;
        if 1==1: #2.9817*(1+1E-2*(x1[0]-1.0)) < 3.1 and 2.9817*(1+1E-2*(x1[0]-1.0)) > 2.9 : #careful. 

            
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
                run.parameters.seed=82*(rn+1)
            #run.parameters.seed=82*(rn+1)    
            #4.02497057364e-05; best critically dressed for 5.2 uT. 103 Hz. 8 significant digits.
            #Best for robust dressing at 5.2 uT 200 Hz mod, B1=4.205e-5; modulation amplitude = 3.497!. 
            #Best for robust dressing at 5.2 uT 300 Hz mod, B1=4.42958e-5; modulation amplitude = 3.32973!. 
            #best for robust dresssing at 5.2 uT 550Hz mod, B1=5.3491388e-05; mod amp= 2.591036

            #best for rcd 5.2 uT 200Hz mod CONSTANT POWER amp= 0.074262070170000, B1=5.22583821E-5
            #8.2968e-05 3225.4323 -2708.6339 0.0011028 0.0036526
            #rcd const power 0.05 opt 8.9866e-05 2897.1981 -2414.8694 0.0012682 0.0034519
            #rcd const power 0.1 opt 8.4395e-05 2966.3542 -2534.4201 0.0012507 0.0034669
            #The reason const power doesn't work is that it isn't floquet in the previous method,
            #The new method optimizes B1, omega, and fm t1 and t2, but solves for wrf_amp so that it is floquet. 
            #hopefully this means it can find a solution at 0.005 and if its good enough it will work for all time. 
            #Hopefully. 
            #at this moment the phase go to 4*Pi for the fm to go one cycle. 
            #it didn't really find a solution. uh ohh. this might not work. 
            #IT DIDNT WORK
            #2, 5.992436587477e+03, 0.0, 1.030520983610e+02, -7.836933948788e+02, 1.084752848683e-03, 6.075986503570e-03, 1.570796327;
            #BzAdd = "3.668263691837e-05";
            #IT FOUND CRITICAL DRESSING!!!!! ugg..
            #not a terrible place to start from?


                wm=2*np.pi*run.parameters.SDparam[3];
                #B1=8.463e-05*(1+1E-2*(x1[0]-1.0));#5.3491388e-05#4.205e-5#5.3491388E-5#4.205e-5#5.347e-05*(1+1E-3*(rn/2-run.parameters.nEntries/4)/(run.parameters.nEntries/2));#4.02444353976e-05*(1+1E-4*(rn/2-run.parameters.nEntries/4)/(run.parameters.nEntries/2));#0.000250021223525*(1+1E-3*(rn/2-run.parameters.nEntries/4)/(run.parameters.nEntries/2))#
            #print "Bz1 = "+str(B1)+"\n"
                #run.EraseFieldFormula();
                #run.parameters.field_assignment("BzAdd",str(B1))
                
                #SDparam[4]=2.9817 is for amplitude RCD
                #run.parameters.SDparam[4]=2.9817*(1+1E-2*(x1[0]-1.0));#3.57E-5/B1*3.497 #*(1+1E-3*(rn/2-run.parameters.nEntries/4)/(run.parameters.nEntries/2));
                #run.parameters.SDparam[1]=2960.3236*(1+1E-2*(x1[1]-1.0));#This is wrf or b  3.57E-5/B1*3.497 #*(1+1E-3*(rn/2-run.parameters.nEntries/4)/(run.parameters.nEntries/2));
                #we fina


                #B1=3.668263691837e-05*(1+1E-2*(x1[0]-1.0));
                #run.EraseFieldFormula();
                #run.parameters.field_assignment("BzAdd",str(B1))

                #run.parameters.SDparam[1]=2.992436587477e+03*(1+1E-2*(x1[1]-1.0));#This is wrf or b  3.57E-5/B1*3.497 #*(1+1E-3*(rn/2-run.parameters.nEntries/4)/(run.parameters.nEntries/2));

                #run.parameters.SDparam[3]=1.030520983610e+02*(1+1E-2*(x1[2]-1.0)); #mo
                #run.parameters.totalTime=1./(run.parameters.SDparam[3]);
                #if run.parameters.SDparam[3] > 600:
                #    run.parameters.SDparam[3] = 600;
                dt0=0.25;#*(1+1E-2*(x1[0]-1.0));
                #dtm=0.017119499912410*(1+1E-2*(x1[0]-1.0));
                dtm=0.040000000000000*(1+1E-2*(x1[0]-1.0));
                #run.parameters.totalTime=2*dt0;
                #run.parameters.SDparam[4]=-2527.5024*(1+1E-2*(x1[2]-1.0));#This is wrf_amp or c 3.57E-5/B1*3.497 #*(1+1E-3*(rn/2-run.parameters.nEntries/4)/(run.parameters.nEntries/2));
                run.parameters.SDparam[5]=dt0-dtm/2.;#0.24*(1+1E-2*(x1[3]-1.0));#this is t1 (transition time) 3.57E-5/B1*3.497 #*(1+1E-3*(rn/2-run.parameters.nEntries/4)/(run.parameters.nEntries/2));
                run.parameters.SDparam[6]=dt0+dtm/2.;#0.26*(1+1E-2*(x1[4]-1.0));#this is t2 (transition time 2) 3.57E-5/B1*3.497 #*(1+1E-3*(rn/2-run.parameters.nEntries/4)/(run.parameters.nEntries/2));
                
                #B1=4.02497057364e-05*(1+1E-2*(x1[1]-1.0));#5.3491388e-05#4.205e-5#5.3491388E-5#4.205e-5#5.347e-05*(1+1E-3*(rn/2-run.parameters.nEntries/4)/(run.parameters.nEntries/2));#4.02444353976e-05*(1+1E-4*(rn/2-run.parameters.nEntries/4)/(run.parameters.nEntries/2));#0.000250021223525*(1+1E-3*(rn/2-run.parameters.nEntries/4)/(run.parameters.nEntries/2))#
            #print "Bz1 = "+str(B1)+"\n"
                #run.EraseFieldFormula();
                #run.parameters.field_assignment("BzAdd",str(B1))
                
                #run.parameters.B0[0]=5.15E-6*(1+1E-2*(x1[0]-1)); 
                

                #if run.parameters.SDparam[5]>0.975*1./(run.parameters.SDparam[3]):
                 #   run.parameters.SDparam[5]=0.975*1./(run.parameters.SDparam[3]);
                  #  print(str(run.parameters.SDparam[5])+"\n")
                #if run.parameters.SDparam[5]<0.025*1./(run.parameters.SDparam[3]):
                 #   run.parameters.SDparam[5]=0.025*1./(run.parameters.SDparam[3]);
                  #  print(str(run.parameters.SDparam[5])+"\n")
            
                #if run.parameters.SDparam[6]>0.975*1./(run.parameters.SDparam[3]):
                 #   run.parameters.SDparam[6]=0.975*1./(run.parameters.SDparam[3]);
                  #  print(str(run.parameters.SDparam[6])+"\n")
                #if run.parameters.SDparam[6]<0.025*1./(run.parameters.SDparam[3]):
                 #   run.parameters.SDparam[6]=0.025*1./(run.parameters.SDparam[3]);
                  #  print(str(run.parameters.SDparam[6])+"\n")
                


            #print "x1 = "+str(run.parameters.SDparam[4])+"\n"
            #print "Mod amp = "+str(run.parameters.SDparam[4]*B1/3.57E-5)+"\n"
            #6 7 Bz1 = 4.1915874e-05

            #
            #run.parameters.Pulseparam[3]=1.245#*(1+.05*(rn/2-run.parameters.nEntries/4)/(run.parameters.nEntries/2))
            #run.parameters.Pulseparam[2]=0.708#*(1+0.1*(rn/2-run.parameters.nEntries/4)/(run.parameters.nEntries/2))#7 is noise rng seed. 
            #Neutrons:
            #3845.309408, 6283.18530718
            #
                phase=np.pi/8.;#0.;#np.pi/2;#np.pi/4.0;
                phi=0.;

                run.parameters.spinSettings[0]=0.
                run.parameters.spinSettings[1]=np.sqrt(1.-np.square(run.parameters.spinSettings[0]))*np.sin(phi+phase)#+phase/2.0)
                run.parameters.spinSettings[2]=np.sqrt(1.-np.square(run.parameters.spinSettings[0]))*np.cos(phi+phase)#+phase/2.0)
                #run.parameters.spinSettings[0]=-0.001441682659922;
                #run.parameters.spinSettings[1]=-0.405022553486012;  
                #run.parameters.spinSettings[2]=0.914305557633103;   


                #run.parameters.edm=1E-23 this is e meters. (ugghhh)
            #run.parameters.edm=0
            #run.parameters.speed=5
            #if rn%4==0:
                #run.parameters.SDparam[2]=np.pi*(rnf/norsf)
                #run.parameters.SDparam[2]=0.0;
                
            #run.parameters.Pulseparam[7]=82*(rn+1) #7 is noise rng seed.
                #print "B1phi = " +str(run.parameters.SDparam[2])+"\n"
                #run.parameters.Pulseparam[1]=(rn//2+2)  
                #offset in another cell, 
                
            else:
                #3He:
                phi=0.;
                phase=np.pi/8.;
                run.parameters.spinSettings[0]=0.
                run.parameters.spinSettings[1]=np.sqrt(1.-np.square(run.parameters.spinSettings[0]))*np.sin(phi-phase)#+phase/2.0)
                run.parameters.spinSettings[2]=np.sqrt(1.-np.square(run.parameters.spinSettings[0]))*np.cos(phi-phase)#+phase/2.0)
                
                #run.parameters.spinSettings[0]=0.001616473758831;
                #run.parameters.spinSettings[1]=0.360218604798280;
                #run.parameters.spinSettings[2]=0.932866519803218;




                #run.parameters.edm=0
                #run.parameters.E
                #print "phi = " +str(run.parameters.SDparam[2])+"\n"
                #7 is noise rng seed. 
                #run.parameters.Pulseparam[7]=82*(rn+1) #7 is noise rng seed. 
             #678742112  #seed is trajectory rng. 
            run.ReloadParameters();
            
            ret = run.runNeutronQuiet();



            if ret == 0:
                ii=0;#run.parameters.nBins-99
                if rn%2==0:
                    #Sxn=run.neutronData.spinx[run.parameters.nBins];
                    #Syn=run.neutronData.spiny[run.parameters.nBins];
                    #Szn=run.neutronData.spinz[run.parameters.nBins];   
                    for i in range(0,run.parameters.nBins+1,1):
                        Sxn[rn//2][ii]=run.neutronData.spinx[i];
                        Syn[rn//2][ii]=run.neutronData.spiny[i];
                        Szn[rn//2][ii]=run.neutronData.spinz[i];
                        tsim[rn//2][ii]=run.neutronData.time[i];
                        #print run.neutronData.spinz[i];
                        ii=ii+1;
                else:
                    #Sx3=run.neutronData.spinx[run.parameters.nBins];
                    #Sy3=run.neutronData.spiny[run.parameters.nBins];
                    #Sz3=run.neutronData.spinz[run.parameters.nBins];
                    for i in range(0,run.parameters.nBins+1, 1):
                        Sx3[rn//2][ii]=run.neutronData.spinx[i];
                        Sy3[rn//2][ii]=run.neutronData.spiny[i];
                        Sz3[rn//2][ii]=run.neutronData.spinz[i];
                        ii=ii+1;
                
                sys.stdout.write(str(rn)+' ')
                sys.stdout.flush()
            
            else:
               print("discarded one neutron")
        else:
            run.nCurrent=run.parameters.nEntries;
    #Done with the calculation, lets do the math to return the value for dumb ol' python.
        
    #print "Sx1 "+str(Sx[0][0])+" Sx2 "+str(Sx[1][0])+"\n";
    #Sxpp=np.mean(Sx,axis=0);
    #print "check if sum of Sx at t = 0?"+str(Sxpp[0])+"\n";
    if 1==1:#(2.9817*(1+1E-2*(x1[0]-1.0)) < 3.1) and (2.9817*(1+1E-2*(x1[0]-1.0)) > 2.9) :    
        #Sxnm=np.mean(Sxn, axis=0);
        #Sznm=np.mean(Szn, axis=0);#Szn[0][1000]##np.mean(Szn, axis=0);
        #Sznm=np.mean(Szn, axis=0);
        #Sx3m=np.mean(Sx3, axis=0);
        #Sz3m=np.mean(Sz3, axis=0);#Sz3[0][1000]#np.mean(Szn, axis=0);#np.mean(Sz3, axis=0);

        #Sxnm=Sxn[0][100];
        #Synm=Syn[0][100];
        #Sznm=Szn[0][100];
        

        #Sx3m=Sx3[0][100];
        #Sy3m=Sy3[0][100];
        #Sz3m=Sz3[0][100];

        tsimm=np.mean(tsim, axis=0);
        #Sz3m=np.mean(Sz3, axis=0);
        Sz32=np.multiply(Sz3,Szn)
        Sy32=np.multiply(Sy3,Syn)              
        Sx32=np.multiply(Sx3,Sxn)
        #S32p=np.add(Sz32,Sy32)
        #S32=np.add(S32p,Sx32)
        phase3=np.arccos(Sy32+Sz32+Sx32);
        phaseMe1=np.mean(phase3[0][0:2000]);
        phaseMe2=np.mean(phase3[0][8000:10000]);
        #Sx=np.square(Sx3m)+np.square(Sxnm);
        #Sz=np.square(Sz3m-1.)+np.square(Sznm-1.0);
        phaseMe=np.square(phaseMe1-phaseMe2);
        #Sznp=np.square(-1.);#np.square(np.mean(Sznm)-1.);
        #Sz3p=np.square(np.mean(Sz3m)-1.);#np.square(np.mean(Sz3m)-1.);
        #print "Sz3= "+str(Sz3p)+", Szn= "+str(Sznp)+"\n";
        #Szn2=np.multiply(Sznm,Sznm)
        #Syn2=np.multiply(Synm,Synm)              
        #Sxn2=np.multiply(Sxnm,Sxnm)
        #Sn2p=np.add(Szn2,Syn2)
        #Sn2=np.add(Sn2p,Sxn2)
        #Sx2=np.abs(Sxn-Sx3);
        #Sy2=np.abs(Syn-Sy3);
        #Sz2=np.abs(Szn-Sz3);
        #S2=np.add(Sznp,Sz3p);
        
        #Sxp=np.subtract(Sxn,Sx3);
        #Syp=np.subtract(Syn,Sy3);
        #Szp=np.subtract(Szn,Sz3);
        #Sx2=np.mean(np.abs(Sxp));
        #Sy2=np.mean(np.abs(Syp));
        #Sz2=np.mean(np.abs(Szp));
        badx=0.0;
        if x1[0]>60.:
            badx=np.square(x1[0]/60.-1)

        S2=phaseMe;#+badx;#Sx+Sz;
        #print("Sxn0 " +str(Sxn[0][0])+" Syn0 "+str(Syn[0][0])+" Syn "+str(Synm)+"  phase_n3 " + str(phasen3)+"\n");
        #print("Sx " +str(Sx)+" Sxn " + str(Sxnm)+" Sx3 " + str(Sx3m)+"\n")# +"#  phase_n3 " + str(phaseMe)+"\n");
        #print "dsx2 = " + str(Sx2)+ " dsy2 = "+str(Sy2)+" dsz2 = "+str(Sz2)+"\n";
        #print("time " + str(tsimm[run.parameters.nBins]) +" x "+str(x1[0])+" "+str(x1[1])+" "+str(x1[2])+" "+str(x1[3])+" "+str(x1[4])+" "+"\n szn = " + str(np.mean(Sznm))+ " sz3 = "+str(np.mean(Sz3m))+"\n");
        #S2 is sigma n cdot sigma 3. THis should be maximized. thus 1-mean(S2) is minimized. 
        #calculating the combined polarization for both neutron and helium3 (this is the maximum signal if a modulation phase is applied).
    
        
        #print str(np.size(S));
    
        #def minExp(x,S):
        #    tfit=np.linspace((run.parameters.nBins-1499)*10.1/2000, 10.1, num=1500)
        #    out1=x[0]*np.exp(-x[1]*tfit);
        #    return np.mean(np.square(np.subtract(out1,S)));



        #res = minimize(minExp, x0, S, method='nelder-mead', options={'xtol': 1e-5, 'disp': True})
        #Sout=2-np.mean(S2);
        print("\n returning value=" +str(S2)+"\n")
        

        return S2;
    else:
        print("\n");
        print("x1 = "+str(x1[0])+" is out of Robust dressing range\n");
        print("returning "+str(np.abs(x1[0]))+"\n");
        return np.abs(x1[0]);








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
    
        
        x0=np.array([1]);#, 1, 1, 1, 1]);

        #x0[0]= 25.7586780717, x0[1]= 1.95572926907, x0[2]= -15.2864911055, x0[3]= -0.904600550417, x0[4]= -0.735150581578
        #x0=np.array([25.7586780717, 1.95572926907,-15.2864911055, -0.904600550417, 0.735150581578]);
       
        #x0=np.array([.95]);
        #def pmin(x1,x2):
            #return (x1-x2)*(x1-x2);
        #x2=1;
        
        def min(x1):
           #del run
           #os.system("rm run"+str(run_id)+".root")
           #run =  create_run(run_id, opts);
           P1m =  single_loop(x1,run);
           return P1m;
        
        #minimize the parameter(s)
        res = minimize(min, x0, method='nelder-mead', options={'xtol': 1e-7, 'disp': True, 'maxiter': 10000, 'maxfev': 10000})
        x0=res.x;

        #print "Optimized x1 =  "+ str(res.x)+"\n"
        print("Saving Optimized Results")
        #print "Running Final Optimization"
        
        #final_loop(x0,run);
        optname = "x1optimized_ConstPow_"+str(run_id);      # optimized value file name
        optd = "/data1/cmswank/spin-sim-xliu/DataSaved/opt/" #optimized x1 value directory
        
        #B1=4.205e-5*(1+1E-2*(x0[1]-1.0));
        
        #x1=2.9817*(1+1E-2*(x0[0]-1.0));


        with open(optd+optname,'w') as foundx1:
            foundx1.writelines(str(x0[0]));#+" , "+str(x0[1]));#+" , "+str(x0[2])+" ,  "+str(x0[3])+" , " +str(x0[4]));
        
        
    
    except KeyboardInterrupt:
        print("^C caught. Letting file finish saving...")
    else:
        print("\ndone")
    finally:
        run.WriteTFile()
        print("exiting..")
        del run

    sys.exit(0)
