import numpy as np 
import pandas as pd 
import rdkit
import cirpy
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw, MolFromSmiles
from rdkit.Chem import PeriodicTable
from rdkit.Chem import GetPeriodicTable
from pkg_resources import resource_filename
import re

names= ['Q','Group','ai','bi','ci','di']
resdata = resource_filename('pyviscosity', 'group_data.csv')
group_data = pd.read_csv(resdata,sep=',', header=0, names=names)
PT = GetPeriodicTable()


class Viscosity(object):

    def __init__(self, inp, T, Pc, debug=False):
        if len(re.findall(r'\b[1-9]{1}[0-9]{1,5}-\d{2}-\d\b',inp)) > 0:
            if debug: print('it\'s a CAS:',inp)
            self.smile = cirpy.resolve(inp,'smiles')
            if debug: print('and it translates to',self.smile)
        else:
            if debug: print('it\'s a smiles:',inp)
            self.smile = inp    
        self.mol= Chem.MolFromSmiles(self.smile)
        print(self.mol)
        self.Qs = []
        self.processed = set()
    
        # here starts the function
        self.atoms = self.mol.GetAtoms()
        self.matrix = rdkit.Chem.rdmolops.GetDistanceMatrix(self.mol)
        self.elements = np.array([PeriodicTable.GetElementSymbol(PT, at.GetAtomicNum()) for at in self.atoms])
    
        for attribute in dir(self):
                if 'check_' in attribute: getattr(self, attribute)()
        #print (elements)
        #print (matrix)
        for Q in self.Qs:
            try: print(group_data[group_data.Q.values==Q][['ai','bi','ci','di']].values[0])
            except: print('None')

    def GetRingSystems(self,includeSpiro=False):
        ri = self.mol.GetRingInfo()
        systems = []
        for ring in ri.AtomRings():
            ringAts = set(ring)
            found = False
            nSystems = []
            for system in systems:
                nInCommon = len(ringAts.intersection(system)) 
                if nInCommon and (includeSpiro or nInCommon>1):
                    ringAts = ringAts.union(system)
                else:
                    nSystems.append(system)
            nSystems.append(ringAts)
            systems = nSystems
        return systems 


    def check_anyhydride(self):
        sub = Chem.MolFromSmarts('C(=O)OC=O')
        match = list(self.mol.GetSubstructMatches(sub))
        for m in match:
            self.processed = self.processed | set(m)
            self.Qs.append(42)
    
    def check_carbonate(self):
        sub = Chem.MolFromSmarts('[O-]C([O-])=O')
        match = list(self.mol.GetSubstructMatches(sub)) 
        for m in match:
            self.processed = self.processed | set(m)
            self.Qs.append(43)

    # #biphenyl 
    # sub = Chem.MolFromSmarts('c1ccc(cc1)c2ccccc2')
    # match = list(mol.GetSubstructMatches(sub)) 
    # for m in match:
    #     processed = processed | set(m)
    #     Qs.append(17) #how many do we add for biphenyl? 
    
    # #terphenyl 
    # sub = Chem.MolFromSmarts('c1ccc(cc1)c2ccc(cc2)c3ccccc3') #what about orpho, meta? 
    # match = list(mol.GetSubstructMatches(sub)) 
    # for m in match:
    #     processed = processed | set(m)
    #     Qs.append(17) #how many do we add for terphenyl?
    
    
    def check_biphenyl_terphenyl_para_meta_ortho(self):
        subs= ['c1ccc(cc1)c2ccccc2','c1ccc(cc1)c2ccc(cc2)c3ccccc3', 'c1ccc(cc1)c2cccc(c2)c3ccccc3', 'c1ccc(cc1)c2ccccc2c3ccccc3']
        for s in subs: 
            sub= Chem.MolFromSmarts(s)
            match = list(self.mol.GetSubstructMatches(sub)) 
            for m in match:
                self.processed = self.processed | set(m)
                self.Qs.append(17) 
    
    #two self.atoms would count as biphenyl, but other atoms in structure should count as e.g. 15/16: need to update logic here 
    
    def check_naphalene(self):
        sub = Chem.MolFromSmarts('c1cccc(c12)cccc2')
        match = list(self.mol.GetSubstructMatches(sub)) 
        for m in match:
            self.processed = self.processed | set(m)
            self.Qs.append(18) #how many do we add for naphalene? 
    
    def check_tetralin(self):
        sub = Chem.MolFromSmarts('C1CCc2ccccc2C1') #this is smile not smarts...
        match = list(self.mol.GetSubstructMatches(sub)) 
        for m in match:
            self.processed = self.processed | set(m)
            self.Qs.append(20) #how many do we add for tetralin? 
    
    def check_aromatic_rings(self):
        sub = Chem.MolFromSmarts("c")
        match = list(self.mol.GetSubstructMatches(sub)) 
        for m in match:
            if m in self.processed : continue
            for el in self.elements[m]:
                indices = np.nonzero(self.matrix[m]==1)[0]
                if self.elements[indices].size == 2: self.Qs.append(15)
                elif self.elements[indices].size == 3: self.Qs.append(16) 
                #add some logic to distinguish between bi/terphenyl, naphthalene, turpentine, tetralin
                #need to know if other carbons in the rings that are not double bonded are classified as this or not 
                #for turpentine: where does classifcation end? leave as cas n.o not recognised anyway 
            self.processed = self.processed | set(m)  
    
    def check_non_aromatic_ring_systems(self):
        rings= self.GetRingSystems()
        for ri in rings: 
            for m in ri: 
                if m in self.processed : continue
                for el in self.elements[m]:
                    hybridization = str(self.atoms[m].GetHybridization())
                    if hybridization == 'SP3': 
                        indices = np.nonzero(self.matrix[m]==1)[0]
                        if self.elements[indices].size == 2:   self.Qs.append(11)
                        elif self.elements[indices].size == 3: self.Qs.append(12)
                        elif self.elements[indices].size == 4: self.Qs.append(14)
                    if hybridization == 'SP2': 
                        if self.elements[indices].size == 1:   self.Qs.append(13)
            self.processed = self.processed | set(ri)  
    
    def check_unlisted(self):
        for i, el in enumerate(self.elements):
            if i in self.processed : continue
            Q = 0
            if el == 'C':
                hybridization = str(self.atoms[i].GetHybridization())
                if hybridization == 'SP3': 
                    indices = np.nonzero(self.matrix[i]==1)[0]
                    if self.elements[indices].size == 0: self.Qs.append(1)
                    elif self.elements[indices].size == 1: self.Qs.append(2)
                    elif self.elements[indices].size == 2: self.Qs.append(3)
                    elif self.elements[indices].size == 3: self.Qs.append(4)
                    elif self.elements[indices].size == 4: self.Qs.append(5)
                elif hybridization == 'SP2': 
                    indices = np.nonzero(self.matrix[i]==1)[0]
                    if self.elements[indices].size == 1: self.Qs.append(6)
                    elif self.elements[indices].size == 2: self.Qs.append(7)
                    elif self.elements[indices].size == 3: self.Qs.append(8)    
                elif hybridization == 'SP': 
                    indices = np.nonzero(self.matrix[i]==1)[0]
                    if self.elements[indices].size == 1: self.Qs.append(9)
                    elif self.elements[indices].size == 2: self.Qs.append(10)
            elif el == 'Cl':
                indices = np.nonzero(self.matrix[i]==2)[0]
                print (self.elements[indices])
                print (self.matrix[i])
                count = np.sum(self.elements[indices]=='Cl')+1
                if   count == 1 : self.Qs.append(78)
                elif count == 2 : self.Qs.append(79)
                elif count == 3 : self.Qs.append(80)
            # TODO processed


#for lab in labels:
 #   print(group_data[group_data.Group.str == lab])


if __name__ == '__main__' :
    import argparse
    parser = argparse.ArgumentParser(
                    prog='pyviscosity',
                    description='What the program does',
                    epilog='Text at the bottom of help')
    parser.add_argument('input')
    print(parser.input)
