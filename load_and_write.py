import myPDB
import collections as col
import h5py
import os
import glob

output = "paired.h5"

files = glob.glob('../PDBPKL4/*_u.pdb.pkl')
matches = col.defaultdict(list)
for file in files:
    pdb_code = os.path.basename(file)[0:4]
    matches[pdb_code].append(file)

for code, files in matches.items():
    if code == '1L9K' or code == 'D35A' or code == 'D35B':
        continue
    if len(files) != 2:
        print "Not enough files for {:}".format(code)
        continue

    files = sorted(files)

    print files
    ligand = myPDB.myPDB.loader(files[0])
    receptor = myPDB.myPDB.loader(files[1])
    with h5py.File(output, 'a') as f:
        grp = f.require_group(code)
        grp['receptor_pssm'] = receptor.pssm.T
        grp['ligand_pssm'] = ligand.pssm.T
        grp['receptor_psfm'] = receptor.psfm.T
        grp['ligand_psfm'] = ligand.psfm.T
        grp['receptor_prasa'] = receptor.rasa
        grp['ligand_prasa'] = ligand.rasa
