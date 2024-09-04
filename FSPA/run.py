# Two files: expected input file + input for multiparameteric
import sys
import os
import subprocess
import shutil


omega0 = 0.05752948549410085 # 0.0568
Ip = 0.5144761644747101 # 0.792482
A01 = 0.1
dA = 0.01
N_pts = 1000

Horders = list(range(13,35,2)) # [37]

program_path = os.environ['FSPA_PATH']+'/FSPA.e'

for order in Horders:
    
    cwd = os.getcwd()
    dname = 'H'+str(order)
    if os.path.exists(dname) and os.path.isdir(dname):
        shutil.rmtree(dname)
        print('deleted previous results')
    os.mkdir(dname)   
    os.chdir(dname)
    content = str(Ip) + '\n' + str(omega0) + '\n' + str(order) + '\n'+ str(A01) + '\n' + \
              str(dA) + '\n' + str(N_pts)
    with open('param.inp','w') as f:
        f.write(content)
        
    subprocess.run(program_path)       
    os.chdir(cwd)