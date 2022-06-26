import requests
import openbabel
def pdb2xyz(pdb_file,type='.pdb'):
    """
    Reads a pdb file and returns the xyz string.
    """
    if type == '.pdb':
        url = "https://files.rcsb.org/view/" + pdb_file + type
    elif type == '.sdf':
        url = "https://files.rcsb.org/ligands/view/" + pdb_file + type 
    else:
        None
    
    if not(None):       
        r = requests.get(url)
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("pdb", "xyz")
        mol = openbabel.OBMol()
        obConversion.ReadString(mol, r.text)
        xyz_mol = obConversion.WriteString(mol)
    return  xyz_mol