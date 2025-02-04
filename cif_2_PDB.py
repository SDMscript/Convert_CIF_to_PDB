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
from biopandas.pdb import PandasPdb

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
       # print (fileName)
        chains =[]
        structure_id = fileName.rsplit('/', 1)[1][:-4]
        joindf=pd.DataFrame()
        mmcif_dict = MMCIF2Dict(fileName)

        entitydf= pd.DataFrame({'id':pd.Series(mmcif_dict['_entity.id']),
                        'entity.description':pd.Series(mmcif_dict['_entity.pdbx_description'])})
        #print (entitydf)
        structdf= pd.DataFrame({'id':pd.Series(mmcif_dict['_struct_ref_seq.ref_id']),
                        'struct_ref_seq.pdbx_strand_id':pd.Series(mmcif_dict['_struct_ref_seq.pdbx_strand_id'])})
        """ _atom_site.group_PDB    ATOM
        _atom_site.id   ATOM_NUM
        _atom_site.type_symbol 
        _atom_site.label_atom_id 
        _atom_site.label_alt_id 
        _atom_site.label_comp_id 
        _atom_site.label_asym_id 
        _atom_site.label_entity_id 
        _atom_site.label_seq_id 
        _atom_site.pdbx_PDB_ins_code 
        _atom_site.Cartn_x  X
        _atom_site.Cartn_y  Y
        _atom_site.Cartn_z  Z
        _atom_site.occupancy    OCC
        _atom_site.B_iso_or_equiv B_ISO
        _atom_site.pdbx_formal_charge 
        _atom_site.auth_seq_id  RES_NUM
        _atom_site.auth_comp_id     RES_NAME
        _atom_site.auth_asym_id     CHAIN
        _atom_site.auth_atom_id     ATOM_NAME
        _atom_site.pdbx_PDB_model_num   """

        coorddata = pd.DataFrame(list(zip(mmcif_dict['_atom_site.group_PDB'], mmcif_dict['_atom_site.id'] , mmcif_dict['_atom_site.type_symbol'], mmcif_dict['_atom_site.Cartn_x'],  mmcif_dict['_atom_site.Cartn_y'],  mmcif_dict['_atom_site.Cartn_z'],  mmcif_dict['_atom_site.occupancy'],  mmcif_dict['_atom_site.B_iso_or_equiv'], mmcif_dict['_atom_site.auth_seq_id'], mmcif_dict['_atom_site.auth_comp_id'], mmcif_dict['_atom_site.label_asym_id'], mmcif_dict['_atom_site.auth_asym_id'], mmcif_dict['_atom_site.auth_atom_id'])),  columns = ['group_PDB', 'id', 'element_id', 'x', 'y', 'z', 'occupancy', 'iso_or_equiv','res_num',  'res_name', 'chain', 'auth_chain','atom_name'])


        for s in segi:
        #    print (s)
            joindf_s=pd.DataFrame()
            t = entitydf[entitydf['entity.description'].str.contains(s)]       
           # print(t)
            joindf_s = pd.merge(structdf,t,on='id')
            joindf_s['segi']=s
            joindf = pd.concat([joindf, joindf_s])
        chains = list(joindf['struct_ref_seq.pdbx_strand_id']) 

        #print(coorddata['auth_asym_id'].unique(), structure_id)
        
        for c in chains:
           # print (structure_id, c)
            segi_chains.append(structure_id+c)
            coordata_chains_c= coorddata[coorddata['auth_chain']   == c]
            coordata_chains_c=coordata_chains_c[coordata_chains_c['group_PDB']=='ATOM']
           # print(coordata_chains_c, structure_id)
            #molidchain =''
            #molidchain = structure_id+c
            seg = joindf[joindf['struct_ref_seq.pdbx_strand_id'] ==c]['segi'].unique()[0]
        #    print (seg)
            coordata_chains_c['segi'] = seg
            #coordata_chains_c['molID_chain'] = molidchain
            coordata_chains_c['molID'] = structure_id
            coorddata_chains = pd.concat([coorddata_chains, coordata_chains_c])
            #print(coordata_chains_c)
            
    return segi_chains, coorddata_chains   
##----
#2. write PDB
def writePDB(coord_dataChains, outname):
    pdbs = list(set(coord_dataChains['molID']))
    for pdb in pdbs:
        coords=coord_dataChains[coord_dataChains['molID']==pdb]
        coords=coords[coords['atom_name'] != 'MG'].reset_index(drop=True)
        #print(coords.head())
        """         # Create a PandasPdb object
        ppdb = PandasPdb()
        # Load the DataFrame into the PandasPdb object
        ppdb.df = coords
        ppdb.to_pdb(path=os.path.join(pdb+".pdb"), 
            records=None, 
            gz=False, 
            append_newline=True) """
        pdb_txt=[]
        for i in range(coords.shape[0]):

            line = "{:6s}{:>5d}  {:3s} {:>3s} {:1s}{:>4d}{:1s}   {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}      {:<3s}\n".format(coords.loc[coords.index[i],'group_PDB'], int(coords.loc[coords.index[i],'id']), coords.loc[coords.index[i],'atom_name'], coords.loc[coords.index[i],'res_name'],  coords.loc[coords.index[i],'chain'],  int(coords.loc[coords.index[i],'res_num']),' ', float(coords.loc[coords.index[i],'x']),  float(coords.loc[coords.index[i],'y']),  float(coords.loc[coords.index[i],'z']),  float(coords.loc[coords.index[i],'occupancy']),  float(coords.loc[coords.index[i],'iso_or_equiv']), coords.loc[coords.index[i],'segi'])
  
            print (line)
            #f.write(line)
            pdb_txt.append(line)
        #    print (line)
        fname = pdb+'.pdb' 
        outfile = os.path.join(outname,fname)
        with open(outfile,'w') as f:
                for l in pdb_txt:
                #    for word in l:
                    f.write(l)
                f.write("TER\n")
                f.write("END\n")

"""     output_pdb = os.path.join(pdb+".pdb")
    pdb_io = PDBIO()
    pdb_io.set_structure(coords)
    pdb_io.save(output_pdb) """



###---------
if __name__ == "__main__":
    ###------
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-f','--folder', help='Path of folder with CIF files to convert to PDB', required=True)
    parser.add_argument('-s','--segi', help='segment keywords to include. For example : \'16S\' for 16S rRNA or 16S ribosomal RNA search. Input can me a comma-separated string (ex. 16S,tRNA).', required=True)
    parser.add_argument('-o','--outfolder', help='Path of folder to save PDB files. Optional. Files will be saved in the source location if output folder is not provided.', required=False)
    
    args = vars(parser.parse_args())
    print(args)
    foldername=args['folder']
    if args['outfolder']:
        outname= args['outfolder']
    else:
        outname=foldername
    segi=list(args['segi'].split(","))
    print (segi)
    pathloc = os.path.abspath(foldername)
    PDB_ssu_chains , coord_dataChains= domain(os.path.abspath(foldername),segi)

    writePDB(coord_dataChains, outname)