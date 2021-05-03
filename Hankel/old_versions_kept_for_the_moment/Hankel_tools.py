### LOAD GRIDS AND FIELDS
print('loading started')

## load radial grid
#file1=open(inpath2+"rgrid.dat","r")
file1=open(os.path.join(inpath2,'rgrid.dat'),"r")
#if file1.mode == "r":
lines = file1.readlines()
k1=0; rgrid=[]
for line in lines: dum = line.split(); rgrid.append(float(dum[0])); k1=k1+1
Nr = k1; rgrid=np.asarray(rgrid);  
file1.close()  

print(rgrid)

## retrieve the dimension of TDSEs
#file1=open( (inpath+"z_000501_r_000001/GridDimensionsForBinaries.dat") ,"r")
file1=open( os.path.join(inpath,'r_000001','GridDimensionsForBinaries.dat') ,"r")
lines = file1.readlines()
file1.close();
Nomega = int(lines[1]);
#print(Nomega)


## omega grid loaded directly in binary form
#file1=open( (inpath+"z_000501_r_000001/omegagrid.bin") ,"rb")
file1=open( os.path.join(inpath,'r_000001','omegagrid.bin') ,"rb")
omegagrid = array.array('d'); #[a.u.]
omegagrid.fromfile(file1,Nomega);
file1.close();

# define corresponding grid points to omegamina and omegamax
#Nomega_anal = 1000 #3000 3000
#Nomega_anal_start = 0
Nomega_anal = FindInterval(omegagrid,omegamax_anal) #3000 3000 1000
Nomega_anal_start = FindInterval(omegagrid,omegamin_anal)
print('om_max', Nomega_anal)
print('om_min', Nomega_anal_start)


## read all fields in the binary form, it follows the padding andf naming of the files
# binary fields are stored in the omega domain (real(1),imaginary(1),real(2),imaginary(2),...)
Nfiles=Nr;
FField_r=np.empty([Nomega,NumOfPointsInRange(0,Nfiles,Nr_step)], dtype=np.cdouble)
k3=0
for k1 in range(0,Nfiles,Nr_step):
  fold=str(k1+1); fold = os.path.join(inpath,'r_'+fold.zfill(6)) # inpath + 'z_000501_r_' + fold.zfill(6) + '/'
  FSourceTermPath = os.path.join(fold,'FSourceTerm.bin')
  file1=open(FSourceTermPath,"rb")
  dum = array.array('d');
  dum.fromfile(file1,2*Nomega);
  for k2 in range(Nomega): FField_r[k2,k3] = dum[2*k2]+1j*dum[2*k2+1];
  k3 = k3+1
#  print(k1)


# reshape rgrid if all points are not used
if ( Nr_step != 1):
  rgridnew=[]
  for k1 in range(int(round(rIntegrationFactorMin*Nfiles)),int(round(rIntegrationFactor*Nfiles)),Nr_step): rgridnew.append(rgrid[k1])
  rgrid=np.asarray(rgridnew);
  Nr = NumOfPointsInRange(int(round(rIntegrationFactorMin*Nfiles)),int(round(rIntegrationFactor*Nfiles)),Nr_step)
  
dr = rgrid[1]-rgrid[0]


print('data loaded')
