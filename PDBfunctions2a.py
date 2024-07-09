import os
import xml.etree.ElementTree as ET
import math
from urllib.request import urlretrieve
import urllib.request
import Bio.PDB
import Bio.PDB.Residue
from alive_progress import alive_it
import shutil
import PDBf_ref_tools as pdbftk
import gzip
import numpy as np
import Bio.Align
from Bio.PDB import MMCIFParser
import Bio.PDB.Chain
import Bio.PDB.Residue
import Bio.PDB.Atom
wdir=os.getcwd()+'/'




class structureFile:

    def __init__(self, path:str,autoseq:bool=False):
        self.filePath=path
        self.title=None
        self.pdbEntryID=None
        self.species=None
        self.chains:dict[str,chain]={}
        #self.ligands={}
        self.hbondCandidates=[]
        self.name=''.join(path.split('/')[-1].split('.')[:-1])
    
        acceptedModes=['cif','mmcif','pdbx','pdb','xml']
        fileExt=path.split('.')[len(path.split('.'))-1]
        fileName=path.split('/')[len(path.split('/'))-1]
        mode=''
        def readPDB():
            file=open(path,'r')
            title=''
            for line in file:
                
                recName=line[0:6].strip()
                if recName == 'ATOM':
                    if len(line) > 80:
                        if len(line) > 82:
                            raise ValueError('"'+fileName+'" cannot be parsed: contains ATOM lines longer than 80 col. wide.')
                        line=line[:28]+line[29:]
                        

                    serial=line[6:11].strip()
                    name=line[12:16].strip()
                    altLoc=line[16].strip()
                    resName=line[17:20].strip()
                    chainID=line[21].strip()
                    resSeq=line[22:26].strip()
                    iCode=line[26].strip()
                    x=line[30:38].strip()
                    y=line[38:46].strip()
                    z=line[46:54].strip()
                    occupancy=line[54:60].strip()
                    tempFactor=line[60:66].strip()
                    element=line[76:78].strip()
                    charge=line[78:80].strip()
                    thisLine=atom(serial,name,altLoc,resName,chainID,resSeq,iCode,x,y,z,occupancy,tempFactor,element,charge)


                    if not chainID in self.chains:
                        self.chains[chainID]=chain(chainID)
                        currentChain=self.chains[chainID]
                        currentChain.residues.append(residue(resSeq, resName))
                        currentResidue=currentChain.residues[-1]

                    if resSeq != currentChain.residues[-1].id:
                        currentChain.residues.append(residue(resSeq, resName))
                        currentResidue=currentChain.residues[-1]
                    #currentChain.atoms.append(thisLine)
                    #currentChain.raw_lines.append(line)
                    currentResidue.atoms[name]=thisLine
                elif recName =='TITLE':
                    title+=line[10:80]
        def readCif():

            parser=MMCIFParser(QUIET=True)

            structure=parser.get_structure(self.filePath.split('/')[-1],self.filePath)

            for model in structure:
                for bio_chain in model:
                    chainIndex=bio_chain.id
                    self.chains[chainIndex] = chain(bio_chain.id,bio_chain)
                    #print(f"Processing chain {bio_chain.id}")

  
        
        def readXML():
            tree=ET.parse(path)
            root=tree.getroot()
            #print (root.attrib)
            tagbase='{http://pdbml.pdb.org/schema/pdbx-v50.xsd}'
            atomLib={
                'label_atom_id':None,       #name
                'label_alt_id':None,        #altLoc
                'label_comp_id':None,       #resName
                'auth_asym_id':None,        #chainID
                'auth_seq_id':None,         #resSeq
                'pdbx_PDB_ins_code':None,   #iCode
                'Cartn_x':None,             #x
                'Cartn_y':None,             #y
                'Cartn_z':None,             #z
                'occupancy':None,           #occupancy
                'B_iso_or_equiv':None,      #tempFactor
                'type_symbol':None,         #element
                'charge':''
                
            }
            for dataBlock in root:
                
                if dataBlock.tag == tagbase+"atom_siteCategory":
                    atoms=dataBlock
                    for entry in atoms:
                        for item in entry:
                            tag=item.tag[len(tagbase):]
                            atomLib[tag]=item.text
                        thisLine=atom(entry.attrib['id'],atomLib['label_atom_id'],atomLib['label_alt_id'],atomLib['label_comp_id'],atomLib['auth_asym_id'],atomLib['auth_seq_id'],atomLib['pdbx_PDB_ins_code'],atomLib['Cartn_x'],atomLib['Cartn_y'],atomLib['Cartn_z'],atomLib['occupancy'],atomLib['B_iso_or_equiv'],atomLib['type_symbol'],atomLib['charge'])

                        if not atomLib['auth_asym_id'] in self.chains:
                            self.chains[atomLib['auth_asym_id']] = chain(atomLib['auth_asym_id'])
                            currentChain=self.chains[atomLib['auth_asym_id']]
                            currentChain.residues.append(residue(atomLib['auth_seq_id'],atomLib['label_comp_id']))
                            currentResidue=currentChain.residues[-1]
                        
                        if atomLib['auth_seq_id'] != currentChain.residues[-1].id:
                            currentChain.residues.append(residue(atomLib['auth_seq_id'],atomLib['label_comp_id']))
                            currentResidue=currentChain.residues[-1]
                        #currentChain.atoms.append(thisLine)
                        #currentChain.raw_lines.append(entry)
                        currentResidue.atoms[atomLib['label_atom_id']]=thisLine
        if fileExt not in acceptedModes:
            raise ValueError('"'+path+'" is not a supported coordinate file type. Supported formats: '+str(acceptedModes).strip('[]'))
        else:
            modeIndex=acceptedModes.index(fileExt)
            if modeIndex <=2:
                mode='cif'
                readCif()
            elif modeIndex == 3:
                mode='legacy'
                readPDB()
            elif modeIndex == 4:
                mode='xml'
                readXML()




class ligand():
    def __init__(self,rawData:Bio.PDB.Residue.Residue):
        self.name:str=rawData.resname
        self.atoms:dict[str,atom] = {}
        pi_sites=[]
        for bio_atom in rawData:
            atom_obj = atom(
            bio_atom.serial_number, bio_atom.name, bio_atom.altloc,
            rawData.resname, rawData.id[0], rawData.id[1],
            rawData.id[2], str(bio_atom.coord[0]), str(bio_atom.coord[1]),
            str(bio_atom.coord[2]), bio_atom.occupancy, bio_atom.bfactor,
            bio_atom.element, str(bio_atom.get_charge())
            )
            self.atoms[bio_atom.id] = atom_obj


class atom():
    def __init__(self,serial:str,name:str,altLoc:str,resName:str,chainID:str,resSeq:str,iCode:str,x:str,y:str,z:str,occupancy:str,tempFactor:str,element:str,charge:str):
        self.serial=int(serial)
        self.name=name
        self.altLoc=altLoc
        self.resName=resName
        self.chainID=chainID
        self.resSeq=int(resSeq)
        self.iCode=iCode
        self.x=float(x)
        self.y=float(y)
        self.z=float(z)
        self.occupancy=float(occupancy)
        self.tempFactor=float(tempFactor)
        self.element=element
        self.charge=charge
    


class residue():
    def __init__(self,rawData:Bio.PDB.Residue.Residue):
        self.id:str=str(rawData.id[1])
        self.type:str=rawData.resname
        self.atoms:dict[str,atom] = {}

        for bio_atom in rawData:
            atom_obj = atom(
            bio_atom.serial_number, bio_atom.name, bio_atom.altloc,
            rawData.resname, rawData.id[0], rawData.id[1],
            rawData.id[2], str(bio_atom.coord[0]), str(bio_atom.coord[1]),
            str(bio_atom.coord[2]), bio_atom.occupancy, bio_atom.bfactor,
            bio_atom.element, str(bio_atom.get_charge())
            )
            self.atoms[bio_atom.id] = atom_obj
            



class chain():
    def __init__(self,name:str,rawData:Bio.PDB.Chain.Chain):
        self.name=name
        self.residues: list[residue] = []
        self.ligands: list[residue] = []
        self.waters: list[residue] = []
        self.__aminosequence=None
        self.distances={}
        self.__bindingSurface=-1.0
        self.__family=None
        self.__type=None
        self.flags=[]
        self.charge_clamps=['h3cc','h4_1cc','h4_8cc','h12cc']
        self.__regPair=None

        
        for bio_residue in rawData:
            if bio_residue.resname in pdbftk.AminoacidDict:
                self.residues.append(residue(bio_residue))
                #print(f"  Adding residue {bio_residue.resname} {bio_residue.id} to chain {self.name}")
            elif bio_residue.resname != 'HOH':
                self.ligands.append(ligand(bio_residue))
                


    
    def seq(self) -> str:
        
        self.__aminosequence = ''
        resBuffer=-1
        for residue in self.residues:
            currentRes=residue.id
            if currentRes != resBuffer:
                self.__aminosequence += pdbftk.AminoacidDict[residue.type]
                resBuffer=currentRes
        
        found=False
        for item in pdbftk.ChrgClmpDict:
            if found: continue
            if pdbftk.ChrgClmpDict[item][0] in self.__aminosequence:
                self.__family='NR'
                self.__type=item
                found=True
        
        if not found:
            for coreg in pdbftk.coacDict:
                if found: continue
                if pdbftk.coacDict[coreg] in self.__aminosequence:
                    self.__family='Coreg'
                    self.__type=coreg
                    found=True

        if not found:
            self.flags.append('unknown')

        return(self.__aminosequence)

    def family(self) -> str:
        if self.__family == None:
            self.seq()
        if 'unknown' in self.flags:
            return 'unknown'
        return self.__family

    def type(self) -> str:
        if self.__type == None:
            self.seq()
        if 'unknown' in self.flags:
            return 'unknown'
        return self.__type
    
    def testForFam(self,test:str) -> bool:
        if self.type() != test:
            return (False)
        return(True)

    
    def cc_dist(self,distReq:list) -> float:
        if self.family() == 'NR':
            if len(distReq) != 2:
                print('Error: only two clamps can be entered to chain.dist()')
                return None
            for clamp in distReq:
                if not clamp in self.charge_clamps:
                    print('Error: "'+clamp+'" is not a valid charge clamp.')
                    print('Valid charge clamps: '+' '.join(self.charge_clamps))
                    return None
            if '|'.join(distReq) not in self.distances:
                

                res1Index = self.seq().index(pdbftk.ChrgClmpDict[self.type()][pdbftk.ccDict[distReq[0]]])
                res2Index = self.seq().index(pdbftk.ChrgClmpDict[self.type()][pdbftk.ccDict[distReq[1]]])

                res1id=self.residues[res1Index].id
                res2id=self.residues[res2Index].id

                distance=self.res_dist([res1id,res2id])
                self.distances['|'.join(distReq)] = distance

            return self.distances['|'.join(distReq)]
        
    def res_dist(self,distReq:list,targetAtom:str='CA'):
        if len(distReq) != 2:
                print('Error: only two clamps can be entered to chain.res_dist()')
                return None
        for res in distReq:
            if type(res) != str:
                raise ValueError('ERROR - requested residues must be the res ID # as type str')
        
        
        if '|'.join(distReq) not in self.distances:
            res1=None
            res2=None
            
            for res in self.residues:
                if res.id == str(distReq[0]):
                    res1=res
                elif res.id == str(distReq[1]):
                    res2=res

            if res1==None or res2==None:
                print('Error: 1 or more residues could not be found')
                return None
            
            res1Atom=res1.atoms[targetAtom]
            res2Atom=res2.atoms[targetAtom]

            distance=atomDist(res1Atom,res2Atom)
            self.distances['|'.join(distReq)] = distance
        
        return self.distances['|'.join(distReq)]

    def bindingSurface(self):
        if self.family() == 'NR':
            if self.__bindingSurface == -1.0:
                reqDistances=['h3cc|h12cc','h3cc|h4_1cc','h4_1cc|h12cc']
                for item in reqDistances:
                    if not item in self.distances:
                        self.cc_dist(item.split('|'))
                
                a=self.distances[reqDistances[0]]
                b=self.distances[reqDistances[1]]
                c=self.distances[reqDistances[2]]
                self.__bindingSurface = pdbftk.heronArea(a,b,c)

            return self.__bindingSurface


def getPDBs(fileList:list, cachePath=wdir+'pdbCache/', cache=True, dlFrmPDB=True) -> dict[str,str]:
    if not cache:
        cachePath = wdir
    if not os.path.exists(cachePath):
        os.mkdir(cachePath)
    
    filePaths={}

    for file in fileList:
        file=file.lower()
        rawFile=file
        if os.path.exists(file):
            filePaths[file] = file
            continue


        if file.endswith('.gz'):
            file=file[:-3]

        if not os.path.exists(cachePath+file):
            if dlFrmPDB and len(file.split('.')[0]) == 4:
                try:
                    url='https://files.rcsb.org/download/'+rawFile
                    urlretrieve(url,cachePath+rawFile)
                    if rawFile.endswith('.gz'):
                        gzipped=gzip.open(cachePath+rawFile,'r')
                        unzipped=open(cachePath+file,'w')
                        gzipcontents=gzipped.readlines()
                        gzipped.close()
                        for line in gzipcontents:
                            unzipped.write(line)
                        unzipped.close()
                except:
                    if os.path.exists(wdir+file):
                        shutil.copyfile(wdir+file,cachePath+file)
                    else:
                        raise ValueError('FILE "'+file+'" DOES NOT EXIST')
        filePaths[file] = cachePath+file

    return filePaths

            


def massInit(codeList:list,extArg='.cif',archivePath=wdir+'pdbCache/',loadingBar=True) -> dict[str,structureFile]:
    filenameList=[]
    for item in codeList:
        if item.strip() == '':
            continue
        ext=extArg
        if item.endswith('.pdb') or item.endswith('.cif') or item.endswith('.xml'):
            ext = ''
        filenameList.append(item+ext)
    filePaths=getPDBs(filenameList, archivePath)
    structures={}
    if loadingBar:
        codeList=alive_it(codeList)
    print('Parsing Structure files...')
    for item in filePaths:
        ext=extArg
        if item.endswith('.pdb') or item.endswith('.cif') or item.endswith('.xml'):
            ext = ''
        iteml=item.lower()
        name=item.split('/')[-1]
        name=''.join(name.split('.')[:-1])
        structures[name]=structureFile(filePaths[iteml+ext],autoseq=True)
    return structures

def atomDist(atom1:atom, atom2:atom) -> float:
    xc1=atom1.x
    yc1=atom1.y
    zc1=atom1.z
    xc2=atom2.x
    yc2=atom2.y
    zc2=atom2.z

    distance=math.sqrt(((round(xc2-xc1,5))**2)+((round(yc2-yc1,5))**2)+((round(zc2-zc1,5))**2))
    return distance

def vectorAngles(p1:list,p2:list,q1:list,q2:list) -> float:
    p_vector=np.array(p2)-np.array(p1)
    q_vector=np.array(q2)-np.array(q1)

    magnitudes=np.linalg.norm(p_vector) * np.linalg.norm(q_vector)
    dp=np.dot(p_vector,q_vector)
    
    cos_theta=np.clip(dp / magnitudes, -1.0, 1.0)

    theta=np.arccos(cos_theta)
    return(np.degrees(theta))


def detectLigandBonding (chain:chain, ligand:ligand, chainRange=-1) -> dict[str,dict[str,any]]:
    
    acceptorAtoms=['N','O','F','S','Se','Cl']
    maxBondDist=3.5
    angleTolerance= 96.0 # In degrees: for all angles which I measured incorrectly + bonding strength
    donors=pdbftk.hBondDonors
    refCarbons=pdbftk.refCarbons
    aproxAngles=pdbftk.aproxAngles
    hbondCandidates={}


    for res in chain.residues:
        for atom in res.atoms:
            if res.atoms[atom].name in donors[pdbftk.AminoacidDict[res.type]]:
                chainatom=res.atoms[atom]
                if not chainatom.name in refCarbons[pdbftk.AminoacidDict[chainatom.resName]]:
                    continue
                for atom1 in ligand.atoms:
                    if ligand.atoms[atom1].element in acceptorAtoms:
                        ligatom=ligand.atoms[atom1]
                        distance=atomDist(chainatom,ligatom)
                        if distance <= maxBondDist:
                            refCarbon=res.atoms[refCarbons[pdbftk.AminoacidDict[chainatom.resName]][chainatom.name]]
                            bondAngle=vectorAngles([refCarbon.x,refCarbon.y,refCarbon.z],[chainatom.x,chainatom.y,chainatom.z],[refCarbon.x,refCarbon.y,refCarbon.z],[ligatom.x,ligatom.y,ligatom.z])
                            idealAngle=aproxAngles[pdbftk.AminoacidDict[chainatom.resName]][chainatom.name]

                            if(bondAngle >= idealAngle-angleTolerance and bondAngle <= idealAngle+angleTolerance):
                                hbondCandidates[str(chainatom.serial)+'|'+str(ligatom.serial)] = {
                                    'residue':chainatom.resSeq,
                                    'atoms':(chainatom,ligatom),
                                    'angle':bondAngle,
                                    'distance':distance
                                }
    
    return hbondCandidates


def findResOffset(chain:chain,MissingResCutoff=math.inf) -> int:
    try:
        upCode=pdbftk.UP_Codes[chain.type()]
    except:
        raise TypeError(f"{chain.name} could not be identified")
    url=f'https://rest.uniprot.org/uniprotkb/{upCode}.fasta'
    fasta=urllib.request.urlopen(url).read().decode('utf-8')
    can_seq=''
    parse=False
    for line in fasta:
        if line=='\n':
            parse=True
            continue
        if line=='>' and parse:
            break
        if parse:
            can_seq+=line.strip()

    #print(can_seq)
    
    #print('Getting Alignment...')
    aligner=Bio.Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.open_gap_score = -1
    aligner.extend_gap_score =-0.1
    aligner.match_score=2
    alignments=aligner.align(chain.seq(),can_seq)

    alignment=alignments[0]

    print(alignment)
    #print(f"score: {alignment.score}")
    #print(alignment[0])
    seqStart=False
    gapOffset=0
    mutations=[]
    for res in range(len(can_seq)):
        align=alignment[0][res]
        if align == '-' and not seqStart:
            continue
        if res < 200:
            continue
        
        if not seqStart:
            firstRes=res
        seqStart = True
        if seqStart:
            if align == '-':
                gapOffset+=1
            elif align == '.':
                mutations.append(res)
        if res >= MissingResCutoff:
            break
    

    index=int(chain.residues[0].id)
    #print(index)
    baseOffset=firstRes-index
    totalOffset= baseOffset+gapOffset
    #print(f'Gap Offset: {gapOffset}')
    #print(f'Base Offset: {baseOffset}')
    print(f'Final Offset: {-totalOffset}')
    return totalOffset