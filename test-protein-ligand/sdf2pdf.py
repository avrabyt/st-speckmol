import requests
from openbabel import openbabel
import streamlit as st

def pdb2xyz(url):
    """
    Reads a pdb file and returns the xyz string.
    """
    
    if not(None):       
        r = requests.get(url)
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("sdf", "pdb")
        mol = openbabel.OBMol()
        obConversion.ReadString(mol, r.text)
        obConversion.WriteFile(mol,"ligand.pdb")
        xyz_mol = obConversion.WriteString(mol)
    return  xyz_mol

pdb_file = st.text_input("Enter pdb file")
if st.button('Convert'):
    st.write(pdb2xyz(pdb_file))
    