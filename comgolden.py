def readaeropara(filename,tokenlist,aeropara):
  mark='----'
  flag=False
  for line in open(filename):
    if(mark in line):
      flag=True
      continue
    if(flag):
      parts=[str for str in line.strip().split(' ') if str!='']
#     print "split line:",parts
      for i in range(len(tokenlist)):
        if(tokenlist[i] in parts[0]):
          if(i<8):
            aeropara.append(float(parts[len(parts)-1]))
          else:
            aeropara.append(float(parts[len(parts)-2]))
            aeropara.append(float(parts[len(parts)-1]))
  return

def buildtoken(tokenlist):
  tokenlist.append('Reynolds')
  tokenlist.append('Mach')
  tokenlist.append('Angle')
  tokenlist.append('Average')
  tokenlist.append('Lift')
  tokenlist.append('Drag')
  tokenlist.append('Friction')
  tokenlist.append('Pitching')
  tokenlist.append('Pressure')
  return
  
import sys

aeropara=[]
aeroparagolden=[]
tokenlist=[]
buildtoken(tokenlist)
#filename='cases/Aeroreport.dat'
filename=sys.argv[1]+'/Aeroreport.dat'
readaeropara(filename,tokenlist,aeropara)
#filename='../2.2.10/ransfoil22-win64-binary/cases/Aeroreport.dat'
#filename='D:\\CaseFiles\\ransfoil\\relaxfactor\\Aeroreport_simple_ra0.7.dat'
filename=sys.argv[2]+'/Aeroreport.dat'
readaeropara(filename,tokenlist,aeroparagolden)
#for a in aeropara:
#  print '%E'%a
#for a in aeroparagolden:
#  print '%E'%a
tol = 1e-2
flag = True
for i in range(len(aeropara)):
  if(aeroparagolden[i] != 0.0):
    res = abs(aeropara[i]-aeroparagolden[i])/aeroparagolden[i]
  else:
    res = abs(aeropara[i]-aeroparagolden[i])
  if(res > tol):
    flag = False
    if(i<9):
      print "The relative error of",tokenlist[i],"with golden is bigger than",tol
    else:
      print "The relative error of",tokenlist[i-1],"with golden is bigger than",tol
if(flag):	  
  print "The relative errors of all the aeropara with golden are smaller than",tol 
