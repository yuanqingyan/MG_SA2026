#!/usr/bin/env python
import os
import sys
import re
import time
import subprocess
import fnmatch
import glob

pCount="./VisiumHD/Apr2025/pThymus_All.txt"
sourcePath="./VisiumHD/Apr2025/Bin2Cell/sourceFile"

tmp='./Spatial/SlurmF'
if not os.path.exists(tmp):
	os.makedirs(tmp)
	
SlurmFolder=tmp
if not os.path.exists(SlurmFolder):
	os.makedirs(SlurmFolder)

def runB2C ():
	sampleName=[]
	slideID=[]
	AreID=[]
	KID=[]
	infile0=open("%s" %pCount,"r").readlines()
	for line in infile0:
		name=line.strip().strip('"').split("\t")[1]
		sampleName.append(name)
		sID=line.strip().strip('"').split("\t")[0]
		slideID.append(sID)
		AID=line.strip().strip('"').split("\t")[2]
		AreID.append(AID)
		keyID=f'{sID}_{AID}_'
		KID.append(keyID)
	print (sampleName)
	
	for i in range(0,len(sampleName)):
		pbsfile=(SlurmFolder+"/rB_%s.slurm" %(sampleName[i]))    
		outfile=open(pbsfile,'w')   
		outfile.write('#!/bin/bash\n')
		outfile.write('#SBATCH -A xxxx\n')     
		outfile.write('#SBATCH -p xxxx\n')   
		outfile.write('#SBATCH -N 1\n')
		outfile.write('#SBATCH --ntasks-per-node=16\n')
		outfile.write('#SBATCH -t 24:00:00\n')
		outfile.write('#SBATCH --mem=180G\n')
		outfile.write('#SBATCH -J B_%s\n'  %(sampleName[i]))
		outfile.write('#SBATCH -o %s/B_%s.o%%j\n' %(SlurmFolder,sampleName[i]))
		outfile.write('#SBATCH -e %s/B_%s.e%%j\n' %(SlurmFolder,sampleName[i]))
		
		outfile.write('#SBATCH --mail-user=xxxx@xxxx\n')
		outfile.write('#SBATCH --mail-type=fail\n')
		outfile.write('date\n')
		outfile.write('module load hdf5/1.8.15-serial\n')
		outfile.write('cd %s\n' %tmp)
		
		
		runB2C=f'python {sourcePath}/S0_sourceBin2Cell.py -ID {sampleName[i]}\n'
		
		
		outfile.write(f"{runB2C}")
		outfile.close()   
	   
		msubCmd=["sbatch %s" %pbsfile]
		subprocess.call(msubCmd,shell=True)   


runB2C()
    


