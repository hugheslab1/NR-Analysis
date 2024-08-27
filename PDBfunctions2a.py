import os
import math
import urllib.request
from alive_progress import alive_it
import shutil
import PDBf_ref_tools as pdbftk
import gzip
import numpy as np
from Bio.PDB import MMCIFParser
from Bio.PDB import PDBParser
import Bio.PDB.Chain
import Bio.PDB.Residue
import Bio.Align
import networkx as nx
import periodictable as ptable
wdir=os.getcwd()+'/'


class structureFile:

    def __init__(self, path:str,autoseq:bool=False):
        self.filePath=path
        self.title=None
        self.pdbEntryID=None
        self.species=None
        self.chains:dict[str,chain]={}
        self.hbondCandidates=[]
        self.name=''.join(path.split('/')[-1].split('.')[:-1])
    
        acceptedModes=['cif','mmcif','pdbx','pdb']
        fileExt=path.split('.')[len(path.split('.'))-1]
        def readPDB():
            parser=PDBParser(QUIET=True)
            structure=parser.get_structure(self.filePath.split('/')[-1],self.filePath)
            for model in structure:
                for bio_chain in model:
                    chainIndex=bio_chain.id
                    self.chains[chainIndex] = chain(bio_chain.id,bio_chain)
        def readCif():
            parser=MMCIFParser(QUIET=True)
            structure=parser.get_structure(self.filePath.split('/')[-1],self.filePath)
            for model in structure:
                for bio_chain in model:
                    chainIndex=bio_chain.id
                    self.chains[chainIndex] = chain(bio_chain.id,bio_chain)

        if fileExt not in acceptedModes:
            raise ValueError(f'"{path}" is not a supported coordinate file type. Supported formats: {str(acceptedModes).strip('[]')}')
        else:
            modeIndex=acceptedModes.index(fileExt)
            if modeIndex <=2:
                mode='cif'
                readCif()
            elif modeIndex == 3:
                mode='legacy'
                readPDB()

        cache = None
        coregPairs={}
        for id in self.chains:
            type = self.chains[id].family()
            
            if cache != None:
                if self.chains[cache].family() != type:
                    coregPairs[cache] = id
                    coregPairs[id] = cache
                elif self.chains[cache].family() == type:
                    cache = id
            elif cache == None:
                cache == id
        self.coregPairs:dict[str,str] = coregPairs
            


class ligand():
    def __init__(self,rawData:Bio.PDB.Residue.Residue):
        self.name:str=rawData.resname
        self.atoms:dict[str,atom] = {}
        self.atomlist:list[atom] = []
        self.__pi_sites:list = []
        self.__ran_pi_scan:bool = False
        for bio_atom in rawData:
            atom_obj = atom(
            bio_atom.serial_number, bio_atom.name, bio_atom.altloc,
            rawData.resname, rawData.id[0], rawData.id[1],
            rawData.id[2], str(bio_atom.coord[0]), str(bio_atom.coord[1]),
            str(bio_atom.coord[2]), bio_atom.occupancy, bio_atom.bfactor,
            bio_atom.element, str(bio_atom.get_charge())
            )
            self.atoms[bio_atom.id] = atom_obj
            self.atomlist.append(atom_obj)

    def __planar(self, ring, tolerance=0.1) -> bool:
        if len(ring) < 3:
            return False
        
        coords = np.array([[self.atoms[atom_name].x, self.atoms[atom_name].y, self.atoms[atom_name].z] for atom_name in ring])

        v1 = coords[1] - coords[0]
        v2 = coords[2] - coords[1]
        
        normal = np.cross(v1,v2)
        normal = normal/np.linalg.norm(normal)

        distances = np.dot(coords-coords[0], normal)

        return np.all(np.abs(distances) < tolerance)

    def pi_sites(self) -> list[list[list]]:
        if not self.__ran_pi_scan:
            G = nx.Graph()
            for atom in self.atomlist:
                G.add_node(atom.name, element=atom.element)
            for i, atom1 in enumerate(self.atomlist):
                for j, atom2 in enumerate(self.atomlist):
                    if i < j and atomDist(atom1, atom2) <= 1.6:
                        G.add_edge(atom1.name, atom2.name)

            cycles=nx.cycle_basis(G)
            for cycle in cycles:
                if self.__planar(cycle) and len(cycle) in [5,6]:
                    #self.__pi_sites.append(cycle)

                    coords = np.array([[self.atoms[atom_name].x, self.atoms[atom_name].y, self.atoms[atom_name].z] for atom_name in cycle])
                    weights = np.array([self.atoms[atom_id].weight for atom_id in cycle])

                    weighted_coords = coords * weights[:, np.newaxis]
                    np_centroid = np.sum(weighted_coords, axis=0) / np.sum(weights)
                    centroid=[]
                    for item in np_centroid:
                        centroid.append(float(item))


                    self.__pi_sites.append([centroid,cycle])
            self.__ran_pi_scan = True
        return self.__pi_sites
    


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
        if len(element) == 2:
            element=element[0]+element[1].lower()
        self.element=element
        self.charge=charge
        self.weight:float = ptable.mass.mass(ptable.elements.isotope(element))
    


class residue():
    def __init__(self,rawData:Bio.PDB.Residue.Residue):
        self.id:str=str(rawData.id[1])
        self.type:str=rawData.resname
        self.atoms:dict[str,atom] = {}
        self.pi_centroids:list[list[float]] = None
        for bio_atom in rawData:
            atom_obj = atom(
            bio_atom.serial_number, bio_atom.name, bio_atom.altloc,
            rawData.resname, rawData.id[0], rawData.id[1],
            rawData.id[2], str(bio_atom.coord[0]), str(bio_atom.coord[1]),
            str(bio_atom.coord[2]), bio_atom.occupancy, bio_atom.bfactor,
            bio_atom.element, str(bio_atom.get_charge())
            )
            self.atoms[bio_atom.id] = atom_obj
        
        if self.type in pdbftk.aa_pi_atoms:
            self.pi_centroids = []
            for loc in pdbftk.aa_pi_atoms[self.type]:
                try: coords = np.array([[self.atoms[atom_name].x, self.atoms[atom_name].y, self.atoms[atom_name].z] for atom_name in loc])
                except: continue
                weights = np.array([self.atoms[atom_id].weight for atom_id in loc])

                weighted_coords = coords * weights[:, np.newaxis]
                np_centroid = np.sum(weighted_coords, axis=0) / np.sum(weights)
                centroid=[]
                for item in np_centroid:
                    centroid.append(float(item))
                
                self.pi_centroids.append(centroid)
            

class chain():
    def __init__(self,name:str,rawData:Bio.PDB.Chain.Chain):
        self.name=name
        self.residues: list[residue] = []
        self.ligands: list[ligand] = []
        self.waters: list[residue] = []
        self.__aminosequence=None
        self.distances={}
        self.__bindingSurface=-1.0
        self.__family=None
        self.__type=None
        self.flags=[]
        self.charge_clamps=['h3cc','h4_1cc','h4_8cc','h12cc']
        self.__regPair=None
        self.__ligand_hydrogen_bonds ={}
        self.__ligand_pi_interactions = {}
        self.__got_h_bonds=False
        self.__got_pi_bonds=False
        
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
                    print(f'Error: "{clamp}" is not a valid charge clamp.')
                    print(f'Valid charge clamps: {" ".join(self.charge_clamps)}')
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
        
    def detectHydrogenBonding (self, maxBondDist=3.5, angleTolerance=96.0) -> dict[str,dict[str,any]]:
        if self.__got_h_bonds == False:
            #acceptorAtoms=['N','O','F','S','Se','Cl']
            nonAcceptors=['C','H']
            maxBondDist=maxBondDist
            angleTolerance= angleTolerance # In degrees: for all angles which I measured incorrectly + bonding strength
            donors=pdbftk.hBondDonors
            refCarbons=pdbftk.refCarbons
            aproxAngles=pdbftk.aproxAngles

            for res in self.residues:
                for atom in res.atoms:
                    if res.atoms[atom].name in donors[pdbftk.AminoacidDict[res.type]]:
                        chainatom=res.atoms[atom]
                        if not chainatom.name in refCarbons[pdbftk.AminoacidDict[chainatom.resName]]:
                            continue
                        for ligand in self.ligands:
                            for atom1 in ligand.atoms:
                                if not ligand.atoms[atom1].element in nonAcceptors:
                                    ligatom=ligand.atoms[atom1]
                                    distance=atomDist(chainatom,ligatom)
                                    if distance <= maxBondDist:
                                        refCarbon=res.atoms[refCarbons[pdbftk.AminoacidDict[chainatom.resName]][chainatom.name]]
                                        bondAngle=vectorAngles([refCarbon.x,refCarbon.y,refCarbon.z],[chainatom.x,chainatom.y,chainatom.z],[refCarbon.x,refCarbon.y,refCarbon.z],[ligatom.x,ligatom.y,ligatom.z])
                                        idealAngle=aproxAngles[pdbftk.AminoacidDict[chainatom.resName]][chainatom.name]

                                        if(bondAngle >= idealAngle-angleTolerance and bondAngle <= idealAngle+angleTolerance):
                                            self.__ligand_hydrogen_bonds[f'{chainatom.serial}|{ligand.name}'] = {
                                                'residue':chainatom.resSeq,
                                                'ligand':ligand.name,
                                                'atoms':(chainatom,ligatom),
                                                'angle':bondAngle,
                                                'distance':distance
                                            }
            self.__got_h_bonds = True

        return self.__ligand_hydrogen_bonds
    
    def detectPiInteractions(self): #Needs more tuning (edge-edge filter, offset edge-side filter, in-between filter)
        if self.__got_pi_bonds == False:
            for res in self.residues:
                if res.type in pdbftk.aa_pi_atoms:
                    for rescentroid in res.pi_centroids:
                        for ligand in self.ligands:
                            for site in ligand.pi_sites():
                                distance = euclidean_distance(rescentroid,site[0])
                                if distance < 6:
                                    #'''
                                    #g=[res.atoms[atom_name].x, res.atoms[atom_name].y, res.atoms[atom_name].z] for atom_name in pdbftk.aa_pi_atoms[res.type][res.pi_centroids.index(rescentroid)]
                                    #print(res.type,res.id)
                                    #print(pdbftk.aa_pi_atoms[res.type][res.pi_centroids.index(rescentroid)])
                                    #print(site[1])
                                    #print([[res.atoms[atom_name].x, res.atoms[atom_name].y, res.atoms[atom_name].z] for atom_name in pdbftk.aa_pi_atoms[res.type][res.pi_centroids.index(rescentroid)]])
                                    #print([[ligand.atoms[atom_name].x, ligand.atoms[atom_name].y, ligand.atoms[atom_name].z] for atom_name in site[1]])

                                    resPlane=fit_plane([[res.atoms[atom_name].x, res.atoms[atom_name].y, res.atoms[atom_name].z] for atom_name in pdbftk.aa_pi_atoms[res.type][res.pi_centroids.index(rescentroid)]])
                                    ligPlane=fit_plane([[ligand.atoms[atom_name].x, ligand.atoms[atom_name].y, ligand.atoms[atom_name].z] for atom_name in site[1]])

                                    #print(resPlane[0])
                                    #print(ligPlane[0])

                                    angle=vector_angle(resPlane[0],ligPlane[0])
                                    #print(angle)
                                    #print()
                                    #print(rescentroid)
                                    #print(angle)
                                    if angle > 90.0:
                                        angle = 180-angle
                                    #print(angle)
                                    go = False

                                    if angle > 30.0:
                                        #print('PING')
                                        res_to_lig = point_to_plane(rescentroid,ligPlane[0],ligPlane[1])
                                        res_to_lig_hyp = float(np.sqrt(res_to_lig**2 + distance**2))
                                        lig_to_res = point_to_plane(site[0],resPlane[0],resPlane[1])
                                        lig_to_res_hyp = float(np.sqrt(lig_to_res**2 + distance**2))
                                        
                                        diffs = [float(np.abs(res_to_lig_hyp - distance)), float(np.abs(lig_to_res_hyp - distance))]
                                        on_plane_ctr  = 0
                                        for item in diffs:
                                            if item > 0.5:
                                                on_plane_ctr += 1
                                        if on_plane_ctr == 1:
                                            go = True
                                       
                                    if angle <= 30.0:
                                        res_to_lig = point_to_plane(rescentroid,ligPlane[0],ligPlane[1])
                                        lig_to_res = point_to_plane(site[0],resPlane[0],resPlane[1])
                                        
                                        if(res_to_lig > 2) and (lig_to_res > 2):
                                            larger=max(res_to_lig,lig_to_res)
                                            smaller=min(res_to_lig,lig_to_res)

                                            if larger * 0.75 < smaller:
                                                go = True

                                    if go:
                                        if angle <= 30.0:
                                            stackType='face-face'
                                        else:
                                            stackType='edge-face'
                                        self.__ligand_pi_interactions[f'{res.id}|{ligand.name}']={
                                            'residue':res.id,
                                            'restype':res.type,
                                            'ligand':ligand.name,
                                            'atoms':site[1],
                                            'distance':distance,
                                            'angle':angle,
                                            'type':stackType
                                        }
                                        
            self.__got_pi_bonds = True
        return self.__ligand_pi_interactions
    

    def ligandBonding(self) -> dict[str,dict[str]]:
        return{
        'hydrogen':self.detectHydrogenBonding(),
        'pi':self.detectPiInteractions()
        }


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
                    url=f'https://files.rcsb.org/download/{rawFile}'
                    urllib.request.urlretrieve(url,cachePath+rawFile)
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
                        raise ValueError(f'FILE "{file}" DOES NOT EXIST')
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

    #print(alignment)
    #print(f"score: {alignment.score}")
    #print(alignment[0])
    seqStart=False
    gapOffset=0
    mutations=[]
    capNoise=0
    for res in range(len(can_seq)):
        align=alignment[0][res]
        if align == '-' and not seqStart:
            continue
        if res < pdbftk.LBD_start[chain.type()]:
            if align != '-':
                capNoise += 1
            continue
        
        if not seqStart:
            firstRes=res + 1
        seqStart = True
        if seqStart:
            if align == '-':
                gapOffset+=1
            elif align == '.':
                mutations.append(res)
        if res >= MissingResCutoff:
            break
    

    index=int(chain.residues[capNoise].id)
    #print(f'index: {index}')
    #print(f'first residue: {firstRes}')
    #print(f'cap noise: {capNoise}')
    baseOffset=firstRes-index
    totalOffset= baseOffset#+gapOffset
    #print(f'Gap Offset: {gapOffset}')
    #print(f'Base Offset: {baseOffset}')
    #print(f'Final Offset: {-totalOffset}')
    return totalOffset

def euclidean_distance(point1:list,point2:list):
    point1 = np.array(point1)
    point2 = np.array(point2)
    distance = np.linalg.norm(point1-point2)
    return float(distance)

def fit_plane(points) -> np.ndarray:
    points=np.array(points)
    A = np.c_[points[:, 0], points[:, 1], np.ones(points.shape[0])]
    b = points[:,2]

    coeffs, _, _, _ = np.linalg.lstsq(A,b,rcond=None)

    normal_vector=[float(coeffs[0]),float(coeffs[1])]

    d = coeffs[2]

    return normal_vector, d

def vector_angle(u:np.ndarray, v:np.ndarray) -> float:


    dot=np.dot(u,v)
    unorm=np.linalg.norm(u)
    vnorm=np.linalg.norm(v)
    
    cos_theta = dot / (unorm * vnorm)
    cos_theta=np.clip(cos_theta,-1.0,1.0)
    
    angle_rad=np.arccos(cos_theta)

    angle_deg=np.degrees(angle_rad)
    return angle_deg

def point_to_plane(point,plane,d) -> float:
    A, B = plane
    x, y, z = point
    distance=(np.abs(A * x + B * y - z + d)) / np.sqrt(A**2 + B**2 + 1)
    #print(distance)
    return(float(distance))

def csvToList(csvPath:str):
    cl=open(csvPath,'r')
    lout=[]
    for line in cl:
        line=line.strip('\n').split(',')
        if len(line) == 1:
            line = line[0]
        lout.append(line)
    cl.close()
    return lout

def checkForExistingOutputFile(filePathName:str,wdir=wdir,rmOld=False):
    i=1
    while os.path.exists(wdir+filePathName):
        filePathName=filePathName.split('(')[0].split('.csv')[0]+'('+str(i)+').csv'
        i+=1
    return(filePathName)


'''
def structures_to_csv(structures:dict[str,structureFile], fileName:str = 'PDBf-out') -> str
    fileName=checkForExistingOutputFile(fileName)    
    out=''
    for structure_id in structures:
        chains = structures[structure_id].chains
        for chain_id in chains:
            chain = chains[chain_id]
            if chain.type() == 'NR':
                out += f'{structures[structure_id].name},{chain.name},'
                if chain_id in structures[structure_id].coregPairs:
                    out += f'{chains[structures[structure_id].coregPairs[chain_id]].type()},'
                else:
                    out += ','
                
                if chain_i
'''