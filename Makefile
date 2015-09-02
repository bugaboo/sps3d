# Start of the makefile Defining variables
objects = sps3d.o pot3d.o potmat.o spsodr.o dvrleg.o dvrr.o legpol.o mul_lrs.o symdvr.o VMULSF.o DLGAMMA.o dvrjac.o indexx.o jacpol.o cbespoln.o beszerr.o kestep.o angbas.o gaujac.o polj.o srad3d.o
FC = gfortran
#f77comp = gfortran
#switch1 = -I/opt/VNI/imsl/fnl600/lnxin100e64/include -L/opt/vni/imsl/fnl600/lnxin100e64/lib -Bdynamic -limsl -limslsuperlu -limslscalar -limslblas -limslmpistub -Xlinker -rpath -Xlinker /opt/vni/imsl/fnl600/lnxin100e64/lib -openmp
switch1 = -L/usr/lib -llapack -lblas -lpthread
switch2 = -c -fdefault-real-8 -fdefault-double-8
# Makefile
sps3d: $(objects)
	gfortran -o sps3d $(objects) $(switch1)
srad3d.o: srad3d.f
	gfortran $(switch2) srad3d.f
sps3d.o: sps3d.f
	gfortran $(switch2) sps3d.f
pot3d.o: pot3d.f
	gfortran $(switch2) pot3d.f
gaujac.o: gaujac.f
	gfortran $(switch2) gaujac.f
polj.o: polj.f
	gfortran $(switch2) polj.f		
potmat.o: potmat.f
	gfortran $(switch2) potmat.f	
spsodr.o: spsodr.f
	gfortran $(switch2) spsodr.f
dvrleg.o: dvrleg.f
	gfortran $(switch2) dvrleg.f
dvrr.o: dvrr.f
	gfortran $(switch2) dvrr.f
legpol.o: legpol.f
	gfortran $(switch2) legpol.f
mul_lrs.o: mul_lrs.f
	gfortran $(switch2) mul_lrs.f
symdvr.o: symdvr.f
	gfortran $(switch2) symdvr.f
VMULSF.o: VMULSF.f
	gfortran $(switch2) VMULSF.f
DLGAMMA.o: DLGAMMA.f
	gfortran $(switch2) DLGAMMA.f
dvrjac.o: dvrjac.f
	gfortran $(switch2) dvrjac.f
indexx.o: indexx.f
	gfortran $(switch2) indexx.f	
jacpol.o: jacpol.f
	gfortran $(switch2) jacpol.f	
#slength.o: slength.f
#	gfortran $(switch2) slength.f	
#spro.o: spro.f
#	gfortran $(switch2) spro.f	
#ssum3d.o: ssum3d.f
#	gfortran $(switch2) ssum3d.f
cbespoln.o: cbespoln.f
	gfortran $(switch2) cbespoln.f
beszerr.o: beszerr.f
	gfortran $(switch2) beszerr.f
#kestep.o: kestep.f
#	gfortran $(switch2) kestep.f
angbas.o : angbas.f
	gfortran $(switch2) angbas.f
	# Cleaning everything
clean:
	rm $(objects)
# End of the makefile
