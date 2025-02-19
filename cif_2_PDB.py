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
    coorddata_chains=pd.DataFrame()
    chains =[]
    structure_id = fileName.rsplit('/', 1)[1][:-4]
    joindf=pd.DataFrame()
    mmcif_dict = MMCIF2Dict(fileName)

    entitydf= pd.DataFrame({'id':pd.Series(mmcif_dict['_entity.id']),
                                'type':pd.Series(mmcif_dict['_entity.type']),
                                'description':pd.Series(mmcif_dict['_entity.pdbx_description']),
                                'pdbx_number_of_molecules':pd.Series(mmcif_dict['_entity.pdbx_number_of_molecules'])})
    if int(entitydf['pdbx_number_of_molecules'][0]) ==1:
        if entitydf.groupby(["id"]).count().max()[0] ==1:
            print (structure_id)
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
            if 'all' in segi:
                    joindf_s=pd.DataFrame()
                    t=entitydf
                    joindf_s = pd.merge(structdf,t,on='id')
                    joindf_s['segi']=t['description'].replace(" ","_")
                    joindf = pd.concat([joindf, joindf_s])

            else:
                    for s in segi:
                #    print (s)
                        joindf_s=pd.DataFrame()
                        t = entitydf[entitydf['description'].str.contains(s, na=False, case=False)]        
                # print(t)
                        joindf_s = pd.merge(structdf,t,on='id')
                        joindf_s['segi']=s
                        joindf = pd.concat([joindf, joindf_s])
                # df = pd.concat([df,joindf])
            if joindf.groupby(["segi"]).count().max()[0] ==1:
                    chains = list(joindf['struct_ref_seq.pdbx_strand_id']) 
                # if joindf.groupby('segi') ==2 :
                    #    print (structure_id, chains)
                    #print(coorddata['auth_asym_id'].unique(), structure_id)
                    
                    for c in chains:
                        print (c)
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
            else:
                    print ('multiple chains '+structure_id)
                    next        #print(coordata_chains_c)        
        else:
                print ('Crystal strucutre contains multiple molecules '+structure_id)
                next        #print(coordata_chains_c)
    else:
            print ('Crystal strucutre contains multiple molecules '+structure_id)
            next
    return segi_chains, coorddata_chains   
##----
#2. write PDB
def writePDB(coord_dataChains, outname):
    pdbs = list(set(coord_dataChains['molID']))
    for pdb in pdbs:
        coords=coord_dataChains[coord_dataChains['molID']==pdb]
        coords=coords[coords['atom_name'] != 'MG'].reset_index(drop=True)
        #print(coords.head())

        pdb_txt=[]
        for i in range(coords.shape[0]):

            line = "{:6s}{:>5d}  {:3s} {:>3s} {:1s}{:>4d}{:1s}   {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}      {:<3s}\n".format(coords.loc[coords.index[i],'group_PDB'], int(coords.loc[coords.index[i],'id']), coords.loc[coords.index[i],'atom_name'], coords.loc[coords.index[i],'res_name'],  coords.loc[coords.index[i],'chain'],  int(coords.loc[coords.index[i],'res_num']),' ', float(coords.loc[coords.index[i],'x']),  float(coords.loc[coords.index[i],'y']),  float(coords.loc[coords.index[i],'z']),  float(coords.loc[coords.index[i],'occupancy']),  float(coords.loc[coords.index[i],'iso_or_equiv']), coords.loc[coords.index[i],'segi'])
  
         #   print (line)
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
    parser.add_argument('-s','--segi', help='default: "all"; segment keywords to include. For example : \'16S\' for 16S rRNA or 16S ribosomal RNA search. Input can me a comma-separated string (ex.  -s "16S,tRNA,mRNA"). all segments will be converted if user does not provide input')
    parser.add_argument('-o','--outfolder', help='Path of folder to save PDB files. Optional. Files will be saved in the source location if output folder is not provided.', required=False)
    
    args = vars(parser.parse_args())
    print(args)
    foldername=args['folder']
    if args['outfolder']:
        outname= args['outfolder']
    else:
        outname=foldername
    if args['segi']:
        segi=list(args['segi'].split(","))
    else:
        segi = 'all'
   # print (segi)
    
    
    pdb_files = glob(os.path.abspath(foldername)+'/*.cif')
    #print(pdb_files)
  #  df = pd.DataFrame()

    for fileName in pdb_files:
        PDB_ssu_chains , coord_dataChains= domain(fileName,segi)
        #df = domain(os.path.abspath(foldername),segi)
        #df.to_csv('./df.csv')
        writePDB(coord_dataChains, outname)