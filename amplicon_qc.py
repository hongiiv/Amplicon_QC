import re
import sys

inFile='/Users/hongiiv/Downloads/amplicon.check/amplicon.info'
data = open(inFile, 'r')
amp={}
amp2={}
for line in data:
   ampchr=line.split('\t')[0]
   ampstart=line.split('\t')[1]
   ampend=line.split('\t')[2]
   ampname=line.split('\t')[3]
   #ampkey='%s_%s_%s'%(ampchr,ampstart,ampend)
   ampkey='%s_%s'%(ampchr,ampstart)
   amp[ampkey]=[ampname]


   amp2chr=line.split('\t')[0]
   amp2start=line.split('\t')[1]
   amp2end=line.split('\t')[2]
   amp2name=line.split('\t')[3]
   #ampkey='%s_%s_%s'%(ampchr,ampstart,ampend)
   amp2key='%s_%s_%s'%(amp2chr,amp2start,amp2end)
   amp2[ampname]=[amp2key]

cigarPattern = '([0-9]+[MIDNSHP])'
cigarSearch = re.compile(cigarPattern)
atomicCigarPattern = '([0-9]+)([MIDNSHP])'
atomicCigarSearch = re.compile(atomicCigarPattern)
softclipPattern = '(^[0-9]+[S])'
softclipSearch = re.compile(softclipPattern)

class queryPos (object):
    """
    struct to store the start and end positions of query CIGAR operations
    """
    def __init__(self, qsPos, qePos, qLen, softLen):
        self.qsPos = int(qsPos)
        self.qePos = int(qePos)
        self.qLen  = int(qLen)
	self.softLen = int(softLen)

class cigarOp (object):
    """
    sturct to store a discrete CIGAR operations
    """
    def __init__(self, opLength, op):
        self.length = int(opLength)
        self.op     = op

def extractCigarOps(cigar):
        if (cigar == "*"):
                cigarOps = []
        else:
                cigarOpStrings = cigarSearch.findall(cigar)
                cigarOps = []
                for opString in cigarOpStrings:
                        cigarOpList = atomicCigarSearch.findall(opString)
                        # "struct" for the op and it's length
                        cigar = cigarOp(cigarOpList[0][0], cigarOpList[0][1])
                        # add to the list of cigarOps
                        cigarOps.append(cigar)
#                       cigarOps = cigarOps
        return(cigarOps)

def calcQueryPosFromCigar(cigarOps):
        qsPos = 0
        qePos = 0
        qLen  = 0
        softLen = 0
        # if first op is a H, need to shift start position
        # the opPosition counter sees if the for loop is looking at the first index of the cigar object
        opPosition = 0
        for cigar in cigarOps:
                if opPosition == 0 and (cigar.op == 'H' or cigar.op == 'S'):
                        qsPos += cigar.length
                        qePos += cigar.length
                        qLen  += cigar.length
			softLen += cigar.length
                elif opPosition > 0 and (cigar.op == 'H' or cigar.op == 'S'):
                        qLen  += cigar.length
                elif cigar.op == 'M' or cigar.op == 'I':
                        qePos += cigar.length
                        qLen  += cigar.length
                        opPosition += 1
        d = queryPos(qsPos, qePos, qLen,softLen);
        return d

qsPos=0
qePos=0
qLen=0
samFile='/Users/hongiiv/Downloads/amplicon.check/HCT15.1.txt'
sam = open(samFile,'r')
for line in sam:
   #print ">>>>"
   samchr=line.split(" ")[2]
   samstart=line.split(" ")[3]
   samcigr=line.split(" ")[5]
   samend=0
   samlen=0

   readCigarOps=extractCigarOps(samcigr)
   readQueryPos = calcQueryPosFromCigar(readCigarOps)
   
   orgsamstart=samstart
   samstart=int(orgsamstart)-readQueryPos.softLen-1
   samend=readQueryPos.qLen+int(orgsamstart)-readQueryPos.softLen-1
   samlen=readQueryPos.qLen

   samkey='%s_%s'%(samchr,samstart)
   samkey_plus_1='%s_%s'%(samchr,int(samstart)+1)
   samkey_minus_1='%s_%s'%(samchr,int(samstart)-1)
   #print samkey
   sys.stdout.write("%s"%line)
   print "Read length: %d"%(readQueryPos.qLen)
   if readQueryPos.softLen > 0:
      print "Clipping length: %d"%readQueryPos.softLen
   print "Read start pos: %s"%(samstart)
   print "Read end pos: %s"%(samend)
   
   matchnum=0
   for j in range(30):
      #print j
      orgsamstart = int(samstart)

      
      samstart_plus = int(samstart)
      samstart_plus = int(samstart_plus)+j
      samkey='%s_%s'%(samchr,samstart_plus)
      
      if amp.get(samkey):
         #print "ooo %d"%j
         amp1=''.join(amp[samkey])
         print "+ Matched Amplicon %d: %s %s"%(matchnum,amp1,amp2[amp1])
         #print "Matched Amplicon: %s"%(amp2[amp1])
         mampstart = str(amp2[amp1]).split("_")[1]
         mampend = str(amp2[amp1]).split("_")[2].strip('\']')
         diffstart = int(mampstart)-int(orgsamstart)
         diffend = int(mampend)-int(samend)
         print "diff start: %d"%(diffstart)
         print "diff end: %d"%(diffend)
         matchnum = matchnum+1

      samstart_minus = int(samstart)-1
      samstart_minus = int(samstart_minus)-j
      samkey='%s_%s'%(samchr,samstart_minus)
      if amp.get(samkey):
         #print "ooo %d"%j
         amp1=''.join(amp[samkey])
         print "- Matched Amplicon %d: %s %s"%(matchnum,amp1,amp2[amp1])
         #print "Matched Amplicon: %s"%(amp2[amp1])
         mampstart = str(amp2[amp1]).split("_")[1]
         mampend = str(amp2[amp1]).split("_")[2].strip('\']')
         diffstart = int(mampstart)-int(orgsamstart)
         diffend = int(mampend)-int(samend)
         print "diff start: %d"%(diffstart)
         print "diff end: %d"%(diffend)
         matchnum = matchnum+1         
   #print "oops %d"%(matchnum) 
   print 

"""
   if amp.get(samkey): 
      #print samkey 
      amp1=''.join(amp[samkey])
      print "Matched Amplicon: %s %s"%(amp1,amp2[amp1])
      #print "Matched Amplicon: %s"%(amp2[amp1])
      mampstart = str(amp2[amp1]).split("_")[1]
      mampend = str(amp2[amp1]).split("_")[2].strip('\']')
      diffstart = int(mampstart)-int(samstart)
      diffend = int(mampend)-int(samend)
      print "diff start: %d"%(diffstart)
      print "diff end: %d"%(diffend)
      #print mampstart
      #print mampend
   elif amp.get(samkey_plus_1):
      #print samkey
      amp1=''.join(amp[samkey_plus_1])
      print "Matched Amplicon: %s %s"%(amp1,amp2[amp1])
      #print "Matched Amplicon: %s"%(amp2[amp1])
      mampstart = str(amp2[amp1]).split("_")[1]
      mampend = str(amp2[amp1]).split("_")[2].strip('\']')
      diffstart = int(mampstart)-int(samstart)
      diffend = int(mampend)-int(samend)
      print "diff start: %d"%(diffstart)
      print "diff end: %d"%(diffend)
   elif amp.get(samkey_minus_1):
      #print samkey
      amp1=''.join(amp[samkey_minus_1])
      print "Matched Amplicon: %s %s"%(amp1,amp2[amp1])
      #print "Matched Amplicon: %s"%(amp2[amp1])
      mampstart = str(amp2[amp1]).split("_")[1]
      mampend = str(amp2[amp1]).split("_")[2].strip('\']')
      diffstart = int(mampstart)-int(samstart)
      diffend = int(mampend)-int(samend)
      print "diff start: %d"%(diffstart)
      print "diff end: %d"%(diffend)
"""

   
