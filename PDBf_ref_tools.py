from math import sqrt

ChrgClmpDict={}
ChrgClmpDict['TRa']=('KKLPMFS','EVFE','EDQII','KGCCM')
ChrgClmpDict['TRb']=('KKLPMFC','EVFE','EDQII','KGCCM')
ChrgClmpDict['RARa']=('KQLPGFT','EMLE','ADQIT','KAACL')
ChrgClmpDict['RARa-worm']=('KQVPVF','EMLD','NDQIT','KAACL')
ChrgClmpDict['RARb']=('KRLPGFTGLT','EMLE','ADQIT','KAACL')
ChrgClmpDict['RARg']=('KRLPGFTGLS','EMLE','ADQIT','KAACL')
ChrgClmpDict['PPARa']=('KAIPGFA','EIYR','NDQVT','KYGVY')
ChrgClmpDict['PPARg']=('KSIPGFV','EIY','NDQVT','KYGVH')
ChrgClmpDict['PPARd']=('KSIPSF','EIYK','NDQVT','KYGVH')
ChrgClmpDict['ReverbAa']=('KHIPGFR','EKLL','HDQVT','KAGTF')
ChrgClmpDict['ReverbAb']=('KRIPGF','EEL','HDQVN','KAGTF')
ChrgClmpDict['RORa']=('KRIDGF','ELF','NDQIV','KAGSL')
ChrgClmpDict['RORb']=('KRITGF','ELFN','NDQIL','KSGCL')
ChrgClmpDict['RORg']=('KRLSGF','ELFS','NDQIV','KAGAM')
ChrgClmpDict['EcR-budworm']=('KGLPGF','EIW','SDQIT','KACSS')
ChrgClmpDict['LXRa-unspecified']=('KQLPGFL','EIWD','EDQIA','KTSAI')
ChrgClmpDict['LXRa-mouse']=('KQLPGFL','EIWD','EDQIA','KTSAI')
ChrgClmpDict['LXRb']=('KQVPGF','EIWD','EDQIA','KASTI')
ChrgClmpDict['FXR']=('KKLPGF','EIW','EDQIA','KGSAV') #HUMAN res99-res269 #Original:EIWDV, exception:Some End after the ASP, not including the Valine;  
ChrgClmpDict['FXR-rat']=('KRLPGFQ','EIWD','EDQIA','KGSAV') #RAT
ChrgClmpDict['VDR-unspecified']=('KMIPGF','EVFG','EDQIV','KSSAI')
ChrgClmpDict['VDR-rat']=('KMIPGF','EVFG','DDQIV','KSSAI')
ChrgClmpDict['VDR-sealamprey']=('KMIPGF','EVFG','EDQIS','KASAI')
ChrgClmpDict['VDR-zebrafish']=('KMIPGF','EVFG','EDQIA','KSSAI')
ChrgClmpDict['PXR']=('KVISYF','ELFG','EDQIS','KGAAF')
ChrgClmpDict['CAR']=('KDLPVF','EICS','EDQIS','KGAAV')
ChrgClmpDict['CAR-mouse']=('KDLPLF','EICS','EDQIS','KGAAV') 
ChrgClmpDict['HNF4a']=('KYIPAF','EML','DDQVA','RAHAG')
ChrgClmpDict['HNF4g']=('KYIPAF','EMLL','DDQVA','RAHAG')
ChrgClmpDict['RXRa']=('KRIPHFSE','EML','DDQVI','RAGWN')
ChrgClmpDict['RXRb']=('KRIPHFSS','EML','DDQVI','RAGWN')
ChrgClmpDict['RXRg']=('KRIPHFSD','EML','EDQVI','RAGWN')
ChrgClmpDict['USP-fruitfly']=('__','__','__','__')
ChrgClmpDict['TR4']=('__','__','__','__')
ChrgClmpDict['TLX-beetle']=('__','__','__','__')
ChrgClmpDict['PNR']=('KNLPVF','LLSD','RDQVI','EEAWS')
ChrgClmpDict['ERa']=('KRVPGF','EMLD','HDQVH','ECAWL')
ChrgClmpDict['ERa-oyster']=('KNVPGY','EMLD','SDQVH','ECCWM') #Missing 'LL' domain after DQ, 'LI' instead
ChrgClmpDict['ERb-unspecified']=('KKIPGF','EML','FDQVR','ESCWM')
ChrgClmpDict['ERb-rat']=('KKIPGF','EMLN','LDQVR','ESCWM')
ChrgClmpDict['ERRa']=('KSIPGFS','EMLE','SDQMS','QSVWM') #Missing 'LL' domain after DQ, 'VL' instead
ChrgClmpDict['ERRb']=('KHIPGFSS','EMLE','GDQMS','QSAWM')
ChrgClmpDict['ERRg-unspecified']=('KHIPGFST','EMLE','ADQMS','QSAWM')
ChrgClmpDict['ERRg-mouse']=('KHIPGFST','EMLE','ADQMS','QSAWM')
ChrgClmpDict['GR-unspecified']=('KAIPGFR','EIIT','DDQMT','QYSWM')
ChrgClmpDict['GR-rat']=('KAIPGFR','EIIT','DDQMT','QYSWM')
ChrgClmpDict['MR']=('KVLPGF','EIIS','EDQIT','QYSWM') #Missing 'LL' domain after DQ, 'LI' instead
ChrgClmpDict['PR']=('KSLPGF','EVIA','DDQIT','QYSWM') #Missing 'LL' domain after DQ, 'LI' instead
ChrgClmpDict['AR-unspecified']=('KALPGF','EIIS','DDQMA','QYSWM') #Missing 'LL' domain after DQ, 'VI' instead
ChrgClmpDict['AR-chimp']=('KALPGF','EIIS','DDQMA','QYSWM')
ChrgClmpDict['AR-rat']=('KALPGF','EIIS','DDQMA','QYSWM')
ChrgClmpDict['AR-mouse']=('KALPGF','EIIS','DDQMA','QYSWM')
ChrgClmpDict['NGFIB']=('EKIPGFAE','KIFM','ADQDL','ESAFL')
ChrgClmpDict['NGFIB-rat']=('EKIPGFI','KIFM','GDQDL','ESAFL')
ChrgClmpDict['NURR1']=('EKIPGFAD','KLFL','ADQDL','ESAFL') #Missing 'LL' domain after DQ, 'LF' instead
ChrgClmpDict['SF1-unspecified']=('__','EMLQ','ADQMT','QNCWS')
ChrgClmpDict['SF1-mouse']=('__','__','__','__')
ChrgClmpDict['LRH-1-unspecified']=('__','__','__','__')
ChrgClmpDict['LRH-1-mouse']=('__','__','__','__')
ChrgClmpDict['Dax-1-mouse']=('__','__','__','__')

helixDict={}
helixDict['FXR']=('ILTEMAT','QVLVEFT','DQIALLK','FLRSAEI')
helixDict['PPARg']=('IRIFQGC','AVQEITE','DQVTLLK','EIIYTML')
helixDict['ERa']=('GLLTNLA','VHMINWA','DQVHLLE','LEILMIG')
helixDict['RXRa']=('NICQAAD','FTLVEWA','DQVILLR','LLIASFS')

coacDict={}
coacDict['SRC1-1']='KLVQLL'
coacDict['SRC1-2']='LHRLLQE'
coacDict['SRC1-3']='QLLRYLL'
coacDict['SRC1-4']='LQQLL'
coacDict['SRC2-1']='LLQLLTT'
coacDict['SRC2-2']='LHRLLQD'
coacDict['SRC2-3']='ALLRYLLDK'
coacDict['SRC3-1']='LLQLLTC'
coacDict['SRC3-2']='ILHKLL'
coacDict['SRC3-3']='LRYLLDR'
coacDict['SRC4-1']='LYSLL'
coacDict['SRC4(b)']='KFKLL'
coacDict['SRC5-1']='LINLL'
coacDict['PGC-1a']='LKKLL'
coacDict['PGC-1b']='LQKLL'
coacDict['NCoR-1']='LGLED'
coacDict['NCoR-1(b)']='LITLAD'
coacDict['NCoR-2']='MGLEA'
coacDict['Med-1']='LMNLL'
coacDict['PA2G4']='LKALL'
coacDict['NRIP1']='LLQLLL'
coacDict['NRIP1(b)']='LLHLL'
coacDict['CREBBP']='LSELL'
coacDict['AR']='FQNLF'
coacDict['SHP']='LKKILL'
coacDict['BUD31']='FDLFY'
coacDict['Q8TBC4']='LQFLL'

ccDict={}
ccDict['h3cc']=0
ccDict['h12cc']=1
ccDict['h4_1cc']=2
ccDict['h4_8cc']=3

AminoacidDict={}
# Alanine
AminoacidDict['ALA']='A' # Alanine
AminoacidDict['DAL']='A' # D-Alanine
AminoacidDict['AKB']='A' # 2-Amino-3-Ketobutyric Acid
AminoacidDict['LPG']='A' # 

AminoacidDict['CYS']='C' # Cysteine
AminoacidDict['CYX']='C' # Cystine (Oxidized Cysteine (PDB - IYY), might change later)
AminoacidDict['CYM']='C' # Cysteine (Deprotonated)
AminoacidDict['ASP']='D' # Aspartic Acid
AminoacidDict['GLU']='E' # Glutamic Acid
AminoacidDict['CGU']='E' # Carboxyglutamic Acid
AminoacidDict['PHE']='F' # Phenylalanine
AminoacidDict['GLY']='G' # Glycine
AminoacidDict['HIS']='H' # Histidine
AminoacidDict['HIE']='H' # Histidine (N3-H tautomer)
AminoacidDict['HID']='H' # Histidine (N1-H tautomer)
AminoacidDict['HIP']='H' # Histidine (Protonated)
AminoacidDict['ILE']='I' # Isoleucine
AminoacidDict['LYS']='K' # Lysine
AminoacidDict['LYZ']='k' # Hydroxylysine
AminoacidDict['ALY']='J' # Acetyllysine **UNOFFICIAL USE OF J**
AminoacidDict['LEU']='L' # Leucine

#Methionine
AminoacidDict['MET']='M' # Methionine
AminoacidDict['MSE']='M' # Selenomethionine

AminoacidDict['ASN']='N' # Asparagine
AminoacidDict['PYL']='O' # Pyrrolysine
AminoacidDict['PRO']='P' # Proline
AminoacidDict['HYP']='P' # Hydroxyproline
AminoacidDict['GLN']='Q' # Glutamine
AminoacidDict['ARG']='R' # Arginine
AminoacidDict['NMM']='r' # Methylarginine/Tilarginine
AminoacidDict['SER']='S' # Serine
AminoacidDict['THR']='T' # Threonine
AminoacidDict['SEC']='U' # Selenocysteine
AminoacidDict['VAL']='V' # Valine
AminoacidDict['TRP']='W' # Tryptophan
AminoacidDict['TYR']='Y' # Tyrosine
AminoacidDict['T44']='Y' # Levothyroxine

#Sometimes Chains are capped with an acetyl or amino group
AminoacidDict['ACE']='@'
AminoacidDict['NH2']='@'
#For DNA in files that contain both DBD and LBD - Just here to avoid key errors
AminoacidDict['DA']='_A_'
AminoacidDict['DT']='_T_'
AminoacidDict['DG']='_G_'
AminoacidDict['DC']='_C_'
AminoacidDict['DU']='_U_'

helixRanges={} # Helix residues for canonical sequences

helixRanges['FXR']={
    1:(260,275),
    2:(280,290),
    3:(293,318),
    4:(326,337),
    5:(338,352),
    6:(354,368),
    7:(372,388),
    8:(392,405),
    9:(414,436),
    10:(442,456),
    11:(457,471),
    12:(476,483)
}



def heronArea (a:float,b:float,c:float):
    s=(a+b+c)/2
    area=sqrt(s*(s-a)*(s-b)*(s-c))
    return area


hBondDonors={
    'A':(),
    'B':(), # OPEN
    'C':('SG'),
    'D':('OD2'),
    'E':('OE2','OE12','OE22'),
    'F':(),
    'G':(),
    'H':('ND1','NE2'),
    'I':(),
    'J':('NZ'), # Unofficial
    'K':('NZ'),
    'k':('OH','NZ'),
    'L':(),
    'M':(),
    'N':('ND2'),
    'O':('NZ'),
    'P':('OD1'),
    'Q':('NE2'),
    'R':('NE','NH2','NH1'),
    'r':('NE','NH1','NH2'),
    'S':('OG'), 
    'T':('OG1'),
    'U':('SE'),
    'V':(),
    'W':('NE1'),
    'X':(),
    'Y':('OH','"O4\'"'),
    'Z':() # OPEN
}

refCarbons={
    'A':{},
    'B':{},
    'C':{
        'SG':'CB'
    },
    'D':{
        'OD2':'CG'
    },
    'E':{
        'OE2':'CD',
        'OE12':'CD1',
        'OE22':'CD2'
    },
    'F':{},
    'G':{},
    'H':{
        'ND1':'CE1',
        'NE2':'CE1'
    },
    'I':{},
    'J':{
        'NZ':'CE'
    },
    'K':{
        'NZ':'CE'
    },
    'k':{
        'OH':'CD',
        'NZ':'CE'
    },
    'L':{},
    'M':{},
    'N':{
        'ND2':'CG'
    },
    'O':{
        'NZ':'CE'
    },
    'P':{
        'OD1':'CG'
    },
    'Q':{
        'NE2':'CD'
    },
    'R':{
        'NE':'CZ',
        'NH2':'CZ',
        'NH1':'CZ',
    },
    'r':{
        'NE':'CZ',
        'NH1':'CZ',
        'NH2':'CZ'
    },
    'S':{
      'OG':'CB'  
    },
    'T':{
        'OG1':'CB'
    },
    'U':{
        'SE':'CB'
    },
    'V':{},
    'W':{
        'NE1':'CD1'
    },
    'X':{},
    'Y':{
        'OH':'CZ',
        '"O4\'"':'"C4\'"'
    }

}

aproxAngles={
    'A':{},
    'B':{},
    'C':{
        'SG':100.0
    },
    'D':{
        'OD2':117.0
    },
    'E':{
        'OE2':117.0,
        'OE12':120.0,
        'OE22':120.0
    },
    'F':{},
    'G':{},
    'H':{
        'ND1':125.0,
        'NE2':125.0
    },
    'I':{},
    'J':{
        'NZ':120.0
    },
    'K':{
        'NZ':109.5
    },
    'k':{
        'OH':107.0,
        'NZ':106.75
    },
    'L':{},
    'M':{},
    'N':{
        'ND2':120.0
    },
    'O':{
        'NZ':120.0
    },
    'P':{
        'OD1':114.0
    },
    'Q':{
        'NE2':120.0
    },
    'R':{
        'NE':120.0,
        'NH2':121.0,
        'NH1':121.0,
    },
    'r':{
        'NE':120.0,
        'NH1':120.0,
        'NH2':120.0
    },
    'S':{
      'OG':107.0  
    },
    'T':{
        'OG1':107.0
    },
    'U':{
        'SE':101.0
    },
    'V':{},
    'W':{
        'NE1':125.0
    },
    'X':{},
    'Y':{
        'OH':107.0,
        '"O4\'"':107.0
    }

}


LigandNames={
    'ADP':'Adenosine Diphosphate',
    'AWL':'XJ-034',
    'AZ2':'AZ-242',
    'BRL':'Rosiglitazone',
    'CHC':'Obeticholic Acid',
    'CTI':'Chelerythrine',
    'C01':'Alpha-Aryloxyphenylacetic Acid Agonist',
    'DKA':'Capric Acid',
    'DRF':'Ragaglitizar',
    'FEZ':'Feroline',
    'GOL':'Glycerol',
    'GWF':'Tropifexor',
    'G24':'GC-24',
    'IVM':'Ivermectin',
    'JN3':'CDCA',
    'KKB':'Seladelpar',
    'KRC':'Cerco-A',
    'LUF':'Luffariellolide',
    'MLA':'Malonic Acid',
    'MLI':'Malonate Ion',
    'MUF':'MFA-1',
    'MYR':'Myristic Acid',
    'NIW':'Nimodipine M (dehydro)',
    'O0U':'BAY-4931',
    'O62':'GSK-8062',
    'REA':'Retinoic Acid',
    'RO0':'PA-082',
    'STI':'Imatinib',
    'TLS':'Telmisartan',
    'TRS':'Tris Buffer',
    'T3':'Liothyronine',
    'T44':'Thyroxine',
    'XD5':'XD-22',
    'XX9':'NDB',
    '0W3':'Ionomycin',
    '034':'GSK-2034',
    '064':'GW-4064',
    '088':'GSK-088',
    '31D':'DM-175',
    '33Y':'XL-335',
    '37G':'GSK-237',
    '486':'Mifepristone',
    '544':'GW-409544',
    '570':'Farglitazar',
    '59G':'GSK-359',
    '8LX':'Lobeglitazone',
    '811':'Cilofexor',
    '82X':'GSK-826',
    '9R0':'HNC-143',
    '9R3':'HNC-180',
    
}

UP_Codes={
    'FXR':'Q96RI1'
}