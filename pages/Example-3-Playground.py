import streamlit as st
from st_speckmol import spec_plot

st.markdown('''# st-speckmol :package:
_A Streamlit **Component** for creating Speck molecular structures within Streamlit Web app._
''')
st.sidebar.header("Add your own xyz coordinates below. :art:")
example_xyz = '''5
methane molecule (in ångströms)
C        0.000000        0.000000        0.000000
H        0.000000        0.000000        1.089000
H        1.026719        0.000000       -0.363000
H       -0.513360       -0.889165       -0.363000
H       -0.513360        0.889165       -0.363000
'''
_xyz = st.sidebar.text_area(
                label = "Paste your coordinates ⬇️",
                value = example_xyz, height  = 200)

st.code(_xyz.splitlines()[1])
res = spec_plot(_xyz)

with st.sidebar.expander("Notes:",expanded=True):
    
    st.markdown(":pushpin:Check out these [xyz_mol_examples folder](https://github.com/avrabyt/st-speckmol/tree/main/xyz_mol_examples). **(The format is critical).**")
    st.markdown("♻️ You can also convert any .pdb files to xyz - [code.](https://github.com/avrabyt/st-speckmol/blob/62d59ff6a059239f64626ef10f0bd9dfa2520d2e/st_speckmol/utils.py#L3)")
    st.markdown(":construction: An [open issue](https://github.com/avrabyt/st-speckmol/issues/7) related to this :arrow_up: feature implementation, feel free to contribute.")