#!/usr/bin/python
import os

outfolder="../experiments/nos5_mid_diffM0"
os.chdir("/home/bgavin/inexactFEAST/bin")      


#lower range for first 15 eigenvalues:
emin=(153670.0+153976.0)/2.0
emax=(161589.0+162059.0)/2
feps=7 #feast tolerance

m0=20
maxloop=200
ellipse=0
int=0
linkspace=1
lineps=16
linits=5
cp=2

ellipses=[100,10000]
m0s=[16,20,30,45]

outputtimes=""
outputfres=""

    #outputtimes=outputtimes+"\n"+str(cp)
    #outputfres=outputfres+"\n"+str(cp)
    
for ellipse in ellipses:
    for m0 in m0s:
        #def runtest(m0,cp,maxloop,ellipse,int,linkspace,linits,lineps):
        ffstr="s\ne !eigenvalue problem type\nd !precision\n"
        ffstr=ffstr+"L !UPLO\n"+str(emin)+" !emin\n"+str(emax)+" !emax\n"
        ffstr=ffstr+str(m0)+" !m0\n11 !#fpm\n1 1 !comments\n2 "+str(cp)+" !contour points\n"
        ffstr=ffstr+"3 "+str(feps)+" !tol\n4 "+str(maxloop)+" !maxloop\n6 1 !convergence\n11 1 !use gmres\n51 "+str(linkspace)
        ffstr=ffstr+" !lin sys krylov size\n50 "+str(linits)+" !gmres iterations \n53 "+str(lineps)+" !lin sys accuracy exponent \n18 "+str(ellipse)+" !ellipse \n16 "+str(int)+" !integration type; 0=gauss, 1=trapezoid, 2=zolotarev"

        #ff = file
        projectname="../matrices/nos5"
        ff=open(projectname+".in","w")
        ff.write(ffstr)
        ff.close()

        basename="nos5_cp"+str(cp)+"_m"+str(m0)+"_int"+str(int)+"_e"+str(ellipse)+"_k"+str(linkspace)+"_r"+str(linits)+"_lineps"+str(lineps)+"_fl"+str(maxloop)+"_feps"+str(feps)

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
