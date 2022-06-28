import streamlit as st
import requests
import openbabel
from st_speckmol import spec_plot


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

# Added from stmol example
st.sidebar.title('Show Proteins')
prot_str='1A2C,1BML,1D5M,1D5X,1D5Z,1D6E,1DEE,1E9F,1FC2,1FCC,1G4U,1GZS,1HE1,1HEZ,1HQR,1HXY,1IBX,1JBU,1JWM,1JWS'
prot_list=prot_str.split(',')
protein=st.sidebar.selectbox('select protein',prot_list)
xyz = pdb2xyz(protein,type='.pdb')
spec_plot(xyz)

st.sidebar.markdown(
    """
    Here the app converts a pdb file to xyz for the speck molecule visualization.
    _The pdb2xyz function is not yet implemented in the st-speckmol library. Will be updated soon._
    """
)
st.sidebar.code('''
xyz = pdb2xyz(protein,type='.pdb'
spec_plot(xyz)
'''
)
