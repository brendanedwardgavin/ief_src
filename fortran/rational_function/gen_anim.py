import os

emin=-0.5
emax=0.5
samples=500
ncp=1

for ellipse in range(0,101):
    os.system("./rational_function "+str(emin)+" "+str(emax)+" "+str(samples)+" "+str(ellipse)+" "+str(ncp))
    ellipse_str=str(ellipse)
    #ellipse_str=ellipse_str.rjust(len(str(100)),'0')
    os.system("cp rational_function.dat output/"+ ellipse_str+".dat")
    os.system("cp contour_points.dat output/"+ ellipse_str+"cp.dat")
    os.system("cp ellipse.dat output/"+ ellipse_str+"el.dat")
