#!/usr/bin/python
import os

def runtest(m0,cp,maxloop,ellipse,int,linkspace,linits,lineps):
    os.chdir("/home/bgavin/inexactFEAST/fortran/sparse_timed")

    ffstr="s\ne !eigenvalue problem type\nd !precision\n"
    ffstr=ffstr+"L !UPLO\n-10.0 !emin\n0.0 !emax\n"
    ffstr=ffstr+str(m0)+" !m0\n11 !#fpm\n1 1 !comments\n2 "+str(cp)+" !contour points\n"
    ffstr=ffstr+"3 10 !tol\n4 "+str(maxloop)+" !maxloop\n6 1 !convergence\n11 1 !use gmres\n51 "+str(linkspace)
    ffstr=ffstr+" !lin sys krylov size\n50 "+str(linits)+" !gmres iterations \n53 "+str(lineps)+" !lin sys accuracy exponent \n18 "+str(ellipse)+" !ellipse \n16 "+str(int)+" !integration type; 0=gauss, 1=trapezoid, 2=zolotarev"

    #ff = file
    projectname="../../matrices/Na5"
    ff=open(projectname+".in","w")
    ff.write(ffstr)
    ff.close()

    basename="na5_cp"+str(cp)+"_m"+str(m0)+"_int"+str(int)+"_e"+str(ellipse)+"_k"+str(linkspace)+"_r"+str(linits)+"_lineps"+str(lineps)

    #basename="na5_m"+str(m0)+"_ss"+str(subspaces)+"_cp"+str(cp)+"_ml"+str(maxloop)+"_e"+str(ellipse)+"_int"+str(int)

    print basename #"m0="+str(m0)+", cp="+str(cp)+", maxloop="+str(maxloop)+", subspaces="+str(subspaces)+", ellipse="+str(ellipse)+", int="+str(int)
    #os.system("./driver_xfeast_sparse "+projectname+" > outfiles/"+basename+".out")
    
    #os.system("../fortran/sparse_timed/time_sddriver "+projectname+" > outfiles/"+basename+".out")
    os.system("./time_sddriver "+projectname+" > ../../autoruns/outfiles/"+basename+".out")

    os.system("cp ../../results/linitsout.dat ../../autoruns/data/"+basename+"_linSysIts.dat")
    os.system("cp ../../results/rhsresidualsout.dat ../../autoruns/data/"+basename+"_linSysRes.dat")
    os.system("cp ../../results/residualsout.dat ../../autoruns/data/"+basename+"_timesEtc.dat")

    ffstr="!!!Input file:\n"+ffstr
    #ffstr=ffstr+"\n\n!!!Output data: iter, contour #, lin_sys parallel, lin_sys seq, residual\n"+feastoutput

    ff=open("../../autoruns/infiles/"+basename+".in","w")
    ff.write(ffstr)
    ff.close()
    return 0


m0=10
cp=3
maxloop=50
ellipse=100
int=0
linkspace=1
linits=100
lineps=2

cp=1
lineps=5
maxloop=150
ellipse=100
runtest(m0,cp,maxloop,ellipse,int,linkspace,linits,lineps)

ellipse=0
runtest(m0,cp,maxloop,ellipse,int,linkspace,linits,lineps)

cp=8
lineps=1
ellipse=100
runtest(m0,cp,maxloop,ellipse,int,linkspace,linits,lineps)

ellipse=0
runtest(m0,cp,maxloop,ellipse,int,linkspace,linits,lineps)













