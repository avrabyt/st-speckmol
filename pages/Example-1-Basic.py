import streamlit as st
import glob
from st_speckmol import speck_plot

st.markdown('''# st-speckmol :package:
_A Streamlit **Component** for creating Speck molecular structures within Streamlit Web app._
''')

# Example files path
ex_files = glob.glob("xyz_mol_examples/*.xyz")
with st.sidebar:
    example_xyz = st.selectbox("Select a molecule",ex_files)
    f = open(example_xyz,"r")
    example_xyz = f.read()
    st.sidebar.info(example_xyz.splitlines()[1])

res = speck_plot(example_xyz)
with st.sidebar.expander('Code', expanded=False):
    st.code(
        '''
        import streamlit as st
import glob
from st_speckmol import speck_plot

# Example files path
ex_files = glob.glob("xyz_mol_examples/*.xyz")
with st.sidebar:
    example_xyz = st.selectbox("Select a molecule",ex_files)
    f = open(example_xyz,"r")
    example_xyz = f.read()
    st.sidebar.info(example_xyz.splitlines()[1])

res = speck_plot(example_xyz)
        '''
    )
