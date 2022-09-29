import mdtraj as md
import streamlit as st

pdb_url = "https://www.rcsb.org/structure/1GFL"
t = md.load_pdb(pdb_url)
st.write(t)


