import streamlit as st
import glob
from st_speckmol import spec_plot

# Example files path
ex_files = glob.glob("xyz_mol_examples/*.xyz")
with st.sidebar:
    example_xyz = st.selectbox("Select a molecule",ex_files)
    f = open(example_xyz,"r")
    example_xyz = f.read()
    st.sidebar.info(example_xyz.splitlines()[1])

res = spec_plot(example_xyz)

