import streamlit as st
import glob
from st_speckmol.utils import *


# Example files path
ex_files = glob.glob("xyz_mol_examples/*.xyz")
with st.sidebar:
    example_xyz = st.selectbox("Select a molecule",ex_files)
    f = open(example_xyz,"r")
    example_xyz = f.read()
    st.sidebar.info(example_xyz.splitlines()[1])

outl = st.sidebar.checkbox('Outline',value=True)
bond = st.sidebar.checkbox('Bond',value=True)
_PARAMETERS = {'outline': outl , 'bondScale': 0.9, 'bonds': bond}
res = spec_plot(example_xyz,_PARAMETERS = _PARAMETERS,wbox_height="500px",wbox_width="500px")
