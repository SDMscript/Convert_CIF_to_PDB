import pandas as pd
from Bio import PDB
from Bio.PDB import *
from glob import glob
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import warnings
import os
import argparse


warnings.filterwarnings('ignore')

# load structures
io = PDBIO()


def domain(fname, segi):
    parser = MMCIFParser()
    segi_chains = []
    coorddata_chains = pd.DataFrame()
    chains = []
    structure_id = fname.rsplit('/', 1)[1][:-4]
    joindf = pd.DataFrame()
    mmcif_dict = MMCIF2Dict(fname)

    entitydf = pd.DataFrame({
        'id': pd.Series(mmcif_dict['_entity.id']),
        'type': pd.Series(mmcif_dict['_entity.type']),
        'description': pd.Series(mmcif_dict['_entity.pdbx_description']),
        'pdbx_number_of_molecules': pd.Series(mmcif_dict['_entity.pdbx_number_of_molecules'])
    })

    if int(entitydf['pdbx_number_of_molecules'][0]) == 1:
        if entitydf.groupby(["id"]).count().max()[0] == 1:
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
                mmcif_dict['_atom_site.label_asym_id'], mmcif_dict['_atom_site.auth_asym_id'],
                mmcif_dict['_atom_site.auth_atom_id']
            )),
                columns=[
                    'group_PDB', 'id', 'element_id', 'x', 'y', 'z', 'occupancy', 'iso_or_equiv',
                    'res_num', 'res_name', 'chain', 'auth_chain', 'atom_name'
                ])

            if 'all' in segi:
                joindf = pd.merge(structdf, entitydf, on='id')
                joindf['segi'] = joindf['struct_ref_seq.pdbx_strand_id']
            else:
                for s in segi:
                    t = entitydf[entitydf['description'].str.contains(s, na=False, case=False)]
                    joindf_s = pd.merge(structdf, t, on='id')
                    joindf_s['segi'] = s
                    joindf = pd.concat([joindf, joindf_s])

            if joindf.groupby(["segi"]).count().max()[0] == 1:
                chains = list(joindf['struct_ref_seq.pdbx_strand_id'])
                for c in chains:
                    segi_chains.append(structure_id + c)
                    coordata_chains_c = coorddata[coorddata['auth_chain'] == c]
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


def writePDB(coord_dataChains, outname):
    pdbs = list(set(coord_dataChains['molID']))
    for pdb in pdbs:
        coords = coord_dataChains[coord_dataChains['molID'] == pdb]
        coords = coords[coords['atom_name'] != 'MG'].reset_index(drop=True)

        pdb_txt = []
        for i in range(coords.shape[0]):
            if len(coords.loc[coords.index[i], 'chain']) == 1:
                line = "{:6s}{:>5d}  {:3s} {:>3s} {:1s}{:>4d}{:1s}   {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}      {:<3s}\n".format(
                    coords.loc[coords.index[i], 'group_PDB'], int(coords.loc[coords.index[i], 'id']),
                    coords.loc[coords.index[i], 'atom_name'], coords.loc[coords.index[i], 'res_name'],
                    coords.loc[coords.index[i], 'chain'], int(coords.loc[coords.index[i], 'res_num']), ' ',
                    float(coords.loc[coords.index[i], 'x']), float(coords.loc[coords.index[i], 'y']),
                    float(coords.loc[coords.index[i], 'z']), float(coords.loc[coords.index[i], 'occupancy']),
                    float(coords.loc[coords.index[i], 'iso_or_equiv']), coords.loc[coords.index[i], 'segi']
                )
            elif len(coords.loc[coords.index[i], 'chain']) == 2:
                next

            pdb_txt.append(line)
        fname = pdb + '.pdb'
        outfile = os.path.join(outname, fname)
        with open(outfile, 'w') as f:
            for l in pdb_txt:
                f.write(l)
            f.write("TER\n")
            f.write("END\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='A Python script to convert multi chain CIF files to PDB format. Allows for segment keyword search. CODE WILL ONLY WRITE CHAINS WITH SINGLE CHARACTER IDENTIFIERS.')
    parser.add_argument('-f', '--folder', help='Path of folder with CIF files to convert to PDB', required=True)
    parser.add_argument('-s', '--segi', help='default: "all"; segment keywords to include. For example : \'16S\' for 16S rRNA or 16S ribosomal RNA search. Input can be a comma-separated string (ex.  -s 16S,rRNA)', required=False)
    parser.add_argument('-o', '--outfolder', help='Path of folder to save PDB files. Optional. Files will be saved in the source location if output folder is not provided.', required=False)

    args = vars(parser.parse_args())
    print(args)
    foldername = args['folder']
    if args['outfolder']:
        outname = args['outfolder']
    else:
        outname = foldername
    if args['segi']:
        segi = list(args['segi'].split(","))
    else:
        segi = 'all'

    pdb_files = glob(os.path.abspath(foldername) + '/*.cif')

    for fileName in pdb_files:
        PDB_seg_chains, coord_dataChains = domain(fileName, segi)
        print(len( PDB_seg_chains), len(coord_dataChains))
        if len( PDB_seg_chains) > 0 and len(coord_dataChains) > 0:
            writePDB(coord_dataChains, outname)
        else:
            print(fileName + ' skipped \n')
