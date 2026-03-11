#!/usr/bin/env python
import os
import sys
import re
import time
import subprocess
import fnmatch
import glob


pCount="./VisiumHD/Apr2025/pThymus_All.txt"
sourcePath="./projects/Manuscript/Thymoma/Fig5__VisiumHD/pythonCode/sourceCode"

tmp='./Spatial/SlurmF'
if not os.path.exists(tmp):
	os.makedirs(tmp)
	
SlurmFolder=tmp
if not os.path.exists(SlurmFolder):
	os.makedirs(SlurmFolder)

def Post1 ():
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
		pbsfile=(SlurmFolder+"/rP1_%s.slurm" %(sampleName[i]))    
		outfile=open(pbsfile,'w')   
		outfile.write('#!/bin/bash\n')
		outfile.write('#SBATCH -A xxx\n')     
		outfile.write('#SBATCH -p xxx\n')   
		outfile.write('#SBATCH -N 1\n')
		outfile.write('#SBATCH --ntasks-per-node=8\n')
		outfile.write('#SBATCH -t 24:00:00\n')
		outfile.write('#SBATCH  --mem=210G\n')
		outfile.write('#SBATCH -J P1_%s\n'  %(sampleName[i]))
		outfile.write('#SBATCH -o %s/P1_%s.o%%j\n' %(SlurmFolder,sampleName[i]))
		outfile.write('#SBATCH -e %s/P1_%s.e%%j\n' %(SlurmFolder,sampleName[i]))
		
		outfile.write('#SBATCH --mail-user=xxx@nxxx\n')
		outfile.write('#SBATCH --mail-type=fail\n')
		
		outfile.write('date\n')
		outfile.write('module load hdf5/1.8.15-serial\n')
		outfile.write('cd %s\n' %tmp)
		
		
		runB2C=f'python {sourcePath}/S3_sourcePostQC.py -ID {sampleName[i]}\n'
		
		
		outfile.write(f"{runB2C}")
		outfile.close()   
	   
		msubCmd=["sbatch %s" %pbsfile]
		subprocess.call(msubCmd,shell=True)   


Post1()
    


