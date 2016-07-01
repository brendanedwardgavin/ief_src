#!/usr/bin/python
import os

projectname="../matrices/parsec/Si2-7e2/Si2"
outfolder="../experiments/Si2_low_diffCP"
os.chdir("/home/bgavin/inexactFEAST/bin")      


#lower range for first 15 eigenvalues:
emin=-1.0
emax=2.0
feps=7 #feast tolerance

m0=23
maxloop=200
ellipse=100
int=0
linkspace=2
lineps=16
linits=5

cps=[1,2,3,4,8,16]

outputtimes=""
outputfres=""

for cp in cps:
    #outputtimes=outputtimes+"\n"+str(cp)
    #outputfres=outputfres+"\n"+str(cp)
    
    
    #def runtest(m0,cp,maxloop,ellipse,int,linkspace,linits,lineps):
    ffstr="s\ne !eigenvalue problem type\nd !precision\n"
    ffstr=ffstr+"L !UPLO\n"+str(emin)+" !emin\n"+str(emax)+" !emax\n"
    ffstr=ffstr+str(m0)+" !m0\n11 !#fpm\n1 1 !comments\n2 "+str(cp)+" !contour points\n"
    ffstr=ffstr+"3 "+str(feps)+" !tol\n4 "+str(maxloop)+" !maxloop\n6 1 !convergence\n11 1 !use gmres\n51 "+str(linkspace)
    ffstr=ffstr+" !lin sys krylov size\n50 "+str(linits)+" !gmres iterations \n53 "+str(lineps)+" !lin sys accuracy exponent \n18 "+str(ellipse)+" !ellipse \n16 "+str(int)+" !integration type; 0=gauss, 1=trapezoid, 2=zolotarev"

    #ff = file

    ff=open(projectname+".in","w")
    ff.write(ffstr)
    ff.close()

    basename="Si2_cp"+str(cp)+"_m"+str(m0)+"_int"+str(int)+"_e"+str(ellipse)+"_k"+str(linkspace)+"_r"+str(linits)+"_lineps"+str(lineps)+"_fl"+str(maxloop)+"_feps"+str(feps)

    #basename="na5_m"+str(m0)+"_ss"+str(subspaces)+"_cp"+str(cp)+"_ml"+str(maxloop)+"_e"+str(ellipse)+"_int"+str(int)

    print basename
    #print "CP = "+str(cp)+", Linits = "+str(linits) #"m0="+str(m0)+", cp="+str(cp)+", maxloop="+str(maxloop)+", subspaces="+str(subspaces)+", ellipse="+str(ellipse)+", int="+str(int)
    #os.system("./driver_xfeast_sparse "+projectname+" > outfiles/"+basename+".out")
    
    #os.system("../fortran/sparse_timed/time_sddriver "+projectname+" > outfiles/"+basename+".out")
    os.system("./time_sddriver "+projectname+" > "+outfolder+"/outfiles/"+basename+".out")

    os.system("cp ../output/linitsout.dat "+outfolder+"/data/"+basename+"_linSysIts.dat")
    os.system("cp ../output/rhsresidualsout.dat "+outfolder+"/data/"+basename+"_linSysRes.dat")
    os.system("cp ../output/residualsout.dat "+outfolder+"/data/"+basename+"_timesEtc.dat")

    ffstr="!!!Input file:\n"+ffstr
    #ffstr=ffstr+"\n\n!!!Output data: iter, contour #, lin_sys parallel, lin_sys seq, residual\n"+feastoutput

    ff=open(outfolder+"/infiles/"+basename+".in","w")
    ff.write(ffstr)
    ff.close()
    
    #finalvalsfile=open("../output/final_vals.dat")
    #filelines=finalvalsfile.readlines()
    #outputtimes=outputtimes+","+filelines[1].strip() #save time
    #outputfres=outputfres+","+filelines[0].strip() #save time
    #finalvalsfile.close()
   
#ff=open(outfolder+"/timesmatrix.dat","w")
#ff.write(outputtimes)
#ff.close()

       
#ff=open(outfolder+"/floopmatrix.dat","w")
#ff.write(outputfres)
#ff.close()
