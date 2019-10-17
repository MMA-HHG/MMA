from scipy import special
import numpy as np
import struct
import array
import os
#import string

def myfun():
  print("abc")


print("Hello world")

myfun()

a = special.jn(0,1)

print(a)

file1=open("35c/fields/rgrid.dat","r")
#if file1.mode == "r":
lines = file1.readlines()
k1=0; z=[]
for xxx in lines: y = xxx.split(); z.append(float(y[1])); k1=k1+1
Nr = k1; z=np.asarray(z);  
file1.close()  


print(z)
#for k1 in range(Nr): print(z[k1])


file1=open("GridDimensionsForBinaries.dat","r")
lines = file1.readlines()
file1.close();
Nomega = int(lines[1]);
print(Nomega)



file1=open("omegagrid.bin","rb")
xxx = array.array('d');
xxx.fromfile(file1,Nomega);
file1.close();
#xxx = struct.unpack('d',file1.read(8))
#print(xxx)

file1=open("FSourceTerm.bin","rb")
xxx = array.array('d');
xxx.fromfile(file1,2*Nomega);
#xxx = struct.unpack('d',file1.read(8))
file1.close();
#print(xxx)
#for k1 in range(Nomega):
#  print(z[k1])

y=np.empty([Nomega], dtype=np.cdouble)
for k1 in range(Nomega): y[k1] = xxx[2*k1]+1j*xxx[2*k1+1];


print(y[Nomega-1])


folders = []

# r=root, d=directories, f = files
#for r, d, f in os.walk("."):
#    for folder in d:
#        folders.append(os.path.join(r, folder))

#for f in folders:
#    print(f)


dirs=os.listdir("35c/TDSEs2")
#print(dirs)
print(dirs[0])


#tstr='xx'
tstr=str(1)
tstr='35c/z_000501_'+tstr.zfill(6)

print(tstr)

Nfiles=2048;
FField_r=np.empty([Nomega,Nfiles], dtype=np.cdouble)
for k1 in range(Nfiles):
  fold=str(k1+1); fold = '35c/z_000501_'+fold.zfill(6)+'/'
  FSourceTermPath = fold+'FSourceTerm.bin'
  file1=open(FSourceTermPath,"rb")
  dum = array.array('d');
  dum.fromfile(file1,2*Nomega);
  for k2 in range(Nomega): FField_r[k2,k1] = dum[2*k2]+1j*dum[2*k2+1];



#print(complex(xxx[1],z[1]))











