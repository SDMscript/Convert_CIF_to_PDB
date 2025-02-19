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

class CIF2PDBConverter:
    def __init__(self, foldername, segi='all', outfolder=None):
        self.foldername = foldername
        self.segi = segi if isinstance(segi, list) else [segi]
        self.outfolder = outfolder if outfolder else foldername
        self.io = PDBIO()
    
    def domain(self, fileName):
        parser = MMCIFParser()
        segi_chains = []
        coorddata_chains = pd.DataFrame()
        chains = []
        structure_id = fileName.rsplit('/', 1)[1][:-4]
        joindf = pd.DataFrame()
        mmcif_dict = MMCIF2Dict(fileName)

        entitydf = pd.DataFrame({
            'id': pd.Series(mmcif_dict['_entity.id']),
            'type': pd.Series(mmcif_dict['_entity.type']),
            'description': pd.Series(mmcif_dict['_entity.pdbx_description']),
            'pdbx_number_of_molecules': pd.Series(mmcif_dict['_entity.pdbx_number_of_molecules'])
        })

        if int(entitydf['pdbx_number_of_molecules'][0]) == 1:
            if entitydf.groupby(["id"]).count().max()[0] == 1:
                print(structure_id)
                structdf = pd.DataFrame({
                    'id': pd.Series(mmcif_dict['_struct_ref_seq.ref_id']),
                    'struct_ref_seq.pdbx_strand_id': pd.Series(mmcif_dict['_struct_ref_seq.pdbx_strand_id'])
                })
                coorddata = pd.DataFrame(list(zip(
                    mmcif_dict['_atom_site.group_PDB'], mmcif_dict['_atom_site.id'],
                    mmcif_dict['_atom_site.type_symbol'], mmcif_dict['_atom_site.Cartn_x'],
                    mmcif_dict['_atom_site.Cartn_y'], mmcif_dict['_atom_site.Cartn_z'],
                    mmcif_dict['_atom_site.occupancy'], mmcif_dict['_atom_site.B_iso_or_equiv'],
                    mmcif_dict['_atom_site.auth_seq_id'], mmcif_dict['_atom_site.auth_comp_id'],
                    mmcif_dict['_atom_site.auth_asym_id'], mmcif_dict['_atom_site.auth_atom_id'],
                    mmcif_dict['_atom_site.pdbx_PDB_model_num']
                )), columns=[
                    'group_PDB', 'id', 'type_symbol', 'Cartn_x', 'Cartn_y', 'Cartn_z',
                    'occupancy', 'B_iso_or_equiv', 'auth_seq_id', 'auth_comp_id',
                    'auth_asym_id', 'auth_atom_id', 'pdbx_PDB_model_num'
                ])
                
                joindf = self._merge_entity_and_struct(entitydf, structdf)

                if joindf.groupby(["segi"]).count().max()[0] == 1:
                    chains = list(joindf['struct_ref_seq.pdbx_strand_id'])
                    for c in chains:
                        print(c)
                        segi_chains.append(structure_id + c)
                        coordata_chains_c = coorddata[coorddata['auth_asym_id'] == c]
                        coordata_chains_c = coordata_chains_c[coordata_chains_c['group_PDB'] == 'ATOM']
                        seg = joindf[joindf['struct_ref_seq.pdbx_strand_id'] == c]['segi'].unique()[0]
                        coordata_chains_c['segi'] = seg
                        coordata_chains_c['molID'] = structure_id
                        coorddata_chains = pd.concat([coorddata_chains, coordata_chains_c])
                else:
                    print('multiple chains ' + structure_id)
                    next
            else:
                print('Crystal structure contains multiple molecules ' + structure_id)
                next
        else:
            print('Crystal structure contains multiple molecules ' + structure_id)
            next
        return segi_chains, coorddata_chains

    def _merge_entity_and_struct(self, entitydf, structdf):
        joindf = pd.DataFrame()
        if 'all' in self.segi:
            joindf_s = pd.DataFrame()
            joindf_s = pd.merge(structdf, entitydf, on='id')
            joindf_s['segi'] = entitydf['description'].replace(" ", "_")
            joindf = pd.concat([joindf, joindf_s])
        else:
            for s in self.segi:
                joindf_s = pd.DataFrame()
                t = entitydf[entitydf['description'].str.contains(s, na=False, case=False)]
                joindf_s = pd.merge(structdf, t, on='id')
                joindf_s['segi'] = s
                joindf = pd.concat([joindf, joindf_s])
        return joindf

    def writePDB(self, coorddata_chains):
        pdbs = list(set(coorddata_chains['molID']))
        for pdb in pdbs:
            coords = coorddata_chains[coorddata_chains['molID'] == pdb]
            coords = coords[coords['auth_atom_id'] != 'MG'].reset_index(drop=True)
            pdb_txt = self._generate_pdb_text(coords)
            self._write_to_file(pdb, pdb_txt)

    def _generate_pdb_text(self, coords):
        pdb_txt = []
        for i in range(coords.shape[0]):
            line = "{:6s}{:>5d}  {:3s} {:>3s} {:1s}{:>4d}{:1s}   {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}      {:<3s}\n".format(
                coords.loc[coords.index[i], 'group_PDB'], int(coords.loc[coords.index[i], 'id']),
                coords.loc[coords.index[i], 'type_symbol'], coords.loc[coords.index[i], 'auth_comp_id'],
                coords.loc[coords.index[i], 'auth_asym_id'], int(coords.loc[coords.index[i], 'auth_seq_id']),
                coords.loc[coords.index[i], 'auth_atom_id'], float(coords.loc[coords.index[i], 'Cartn_x']),
                float(coords.loc[coords.index[i], 'Cartn_y']), float(coords.loc[coords.index[i], 'Cartn_z']),
                float(coords.loc[coords.index[i], 'occupancy']), float(coords.loc[coords.index[i], 'B_iso_or_equiv']),
                coords.loc[coords.index[i], 'pdbx_PDB_model_num']
            )
            pdb_txt.append(line)
        return pdb_txt

    def _write_to_file(self, pdb, pdb_txt):
        fname = pdb + '.pdb'
        outfile = os.path.join(self.outfolder, fname)
        with open(outfile, 'w') as f:
            for l in pdb_txt:
                f.write(l)
            f.write("TER\n")
            f.write("END\n")

    def convert(self):
        pdb_files = glob(os.path.abspath(self.foldername) + '/*.cif')
        for fileName in pdb_files:
            PDB_ssu_chains, coord_dataChains = self.domain(fileName)
            self.writePDB(coord_dataChains)
            
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-f', '--folder', help='Path of folder with CIF files to convert to PDB', required=True)
    parser.add_argument('-s', '--segi', help='default: "all"; segment keywords to include. For example: "16S" for 16S rRNA or 16S ribosomal RNA search. Input can be a comma-separated string (ex. -s 16S,rRNA)', required=False, default='all')
    parser.add_argument('-o', '--outfolder', help='Path of folder to save PDB files. Optional. Files will be saved in the source location if output folder is not provided.', required=False)
    
    args = vars(parser.parse_args())
    foldername = args['folder']
    outname = args['outfolder'] if args['outfolder'] else foldername
    segi = list(args['segi'].split(",")) if args['segi'] else 'all'
    
    converter = CIF2PDBConverter(foldername, segi, outname)
    converter.convert()