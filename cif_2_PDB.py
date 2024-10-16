import pandas as pd
from Bio import PDB
from Bio.PDB import *
from glob import glob
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import warnings
warnings.filterwarnings('ignore')
import os
from os.path import join
import argparse

#test 

###------FUNCTION--------
#### 1. LOAD CIF AND EXTRACT DATA

# #load structures
io = PDBIO()

def domain(fname,segi):
    parser = MMCIFParser()
    segi_chains =[]
    pdb_files = glob(fname+'/*.cif')
    #print(pdb_files)
    coorddata_chains=pd.DataFrame()
    for fileName in pdb_files:
        #print (fileName)
        chains =[]
        structure_id = fileName.rsplit('/', 1)[1][:-4]
        joindf=pd.DataFrame()
        mmcif_dict = MMCIF2Dict(fileName)

        entitydf= pd.DataFrame({'id':pd.Series(mmcif_dict['_entity.id']),
                        'entity.description':pd.Series(mmcif_dict['_entity.pdbx_description'])})
        #print (entitydf)
        structdf= pd.DataFrame({'id':pd.Series(mmcif_dict['_struct_ref_seq.ref_id']),
                        'struct_ref_seq.pdbx_strand_id':pd.Series(mmcif_dict['_struct_ref_seq.pdbx_strand_id'])})
        
        coorddata = pd.DataFrame(list(zip(mmcif_dict['_atom_site.group_PDB'], mmcif_dict['_atom_site.id'],  mmcif_dict['_atom_site.type_symbol'],  mmcif_dict['_atom_site.label_atom_id'],  mmcif_dict['_atom_site.label_alt_id'], 
                                        mmcif_dict['_atom_site.label_comp_id'],  mmcif_dict['_atom_site.auth_asym_id'],  mmcif_dict['_atom_site.label_entity_id'],  mmcif_dict['_atom_site.auth_seq_id'],  mmcif_dict['_atom_site.pdbx_PDB_ins_code'], 
                                        mmcif_dict['_atom_site.Cartn_x'],  mmcif_dict['_atom_site.Cartn_y'],  mmcif_dict['_atom_site.Cartn_z'],  mmcif_dict['_atom_site.occupancy'],  mmcif_dict['_atom_site.B_iso_or_equiv'])),
                                        columns = ['group_PDB', 'id', 'type_symbol', 'label_atom_id', 'label_alt_id', 
                                        'label_comp_id', 'auth_asym_id', 'label_entity_id', 'auth_seq_id', 'pdbx_PDB_ins_code', 
                                        'Cartn_x', 'Cartn_y', 'Cartn_z','occupancy',  'B_iso_or_equiv'])
        for s in segi:
            joindf_s=pd.DataFrame()
            t = entitydf[entitydf['entity.description'].str.contains(s)]       
            print(t)
            joindf_s = pd.merge(structdf,t,on='id')
            joindf_s['segi']=s
            joindf = pd.concat([joindf, joindf_s])
        chains = list(joindf['struct_ref_seq.pdbx_strand_id']) 

        #print(coorddata['auth_asym_id'].unique(), structure_id)
        
        for c in chains:
           # print (structure_id, c)
            segi_chains.append(structure_id+c)
            coordata_chains_c= coorddata[coorddata['auth_asym_id']   == c]
           # print(coordata_chains_c, structure_id)
            #molidchain =''
            #molidchain = structure_id+c
            seg = joindf[joindf['struct_ref_seq.pdbx_strand_id'] ==c]['segi'].unique()[0]
            print (seg)
            coordata_chains_c['segi'] = seg
            #coordata_chains_c['molID_chain'] = molidchain
            coordata_chains_c['molID'] = structure_id
            coorddata_chains = pd.concat([coorddata_chains, coordata_chains_c])
            #print(coordata_chains_c)
            
    return segi_chains, coorddata_chains   
##----
#2. write PDB
def writePDB(coord_dataChains):
    pdbs = list(set(coord_dataChains['molID']))
    for pdb in pdbs:
        coords=coord_dataChains[coord_dataChains['molID']==pdb]
        print(coords.head())
        
        pdb_txt=[]
        for i in range(coords.shape[0]):
            #print([coords.iloc[i,2]][0])
        # print (coords.iloc[i,18])
            #coords.iloc[i,17]
            line = "{:6s}{:5d}{:>4s}{:>5s} {:>1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}      {:3s}".format(coords.iloc[i,0],int(coords.iloc[i,1]),coords.iloc[i,3],coords.iloc[i,5],coords.iloc[i,6],int(coords.iloc[i,8]),' ',float(coords.iloc[i,10]),float(coords.iloc[i,11]),float(coords.iloc[i,12]), float(coords.iloc[i,13]),float(coords.iloc[i,14]),coords.iloc[i,15])
            pdb_txt.append(line+"\n")
            #print (line)
        fname = pdb+'_RNA.pdb' 
        outfile = os.path.join(foldername,fname)
        with open(outfile,'w') as f:
                for l in pdb_txt:
                #    for word in l:
                    f.write(l)
                f.write("TER\n")
                f.write("END\n")


###---------
if __name__ == "__main__":
    ###------
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-f','--folder', help='Path of folder with CIF files to convert to PDB', required=True)
    parser.add_argument('-s','--segi', help='segment keywords to include. For example : \'16S\' for 16S rRNA or 16S ribosomal RNA search. Input can me a comma-separated list or string.', required=True)
    args = vars(parser.parse_args())
    print(args)
    foldername=args['folder']
    segi=args['segi']
    pathloc = os.path.abspath(foldername)
    PDB_ssu_chains , coord_dataChains= domain(os.path.abspath(foldername),segi)

    writePDB(coord_dataChains)