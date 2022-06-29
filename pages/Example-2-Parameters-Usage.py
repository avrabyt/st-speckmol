import streamlit as st
import glob
from st_speckmol import spec_plot

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

with st.sidebar.expander("Parameters",expanded=True):
    outl = st.checkbox('Outline',value=True)
    bond = st.checkbox('Bond',value=True)
    bond_scale = st.slider('BondScale',min_value=0.0,max_value=1.0,value=0.8)
    brightness = st.slider('Brightness',min_value=0.0,max_value=1.0,value=0.4)
    relativeAtomScale = st.slider('RelativeAtomScale',min_value=0.0,max_value=1.0,value=0.64)
    bondShade = st.slider('bondShade',min_value=0.0,max_value=1.0,value=0.5)

_PARAMETERS = {'outline': outl , 'bondScale': bond_scale,
                'bonds': bond ,'bondShade' : bondShade,
                'brightness': brightness, 'relativeAtomScale':relativeAtomScale,
                }
res = spec_plot(example_xyz,_PARAMETERS = _PARAMETERS,wbox_height="500px",wbox_width="500px")

with st.expander('Code',expanded=False):
    st.code(
        '''
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

    with st.sidebar.expander("Parameters",expanded=True):
        outl = st.checkbox('Outline',value=True)
        bond = st.checkbox('Bond',value=True)
        bond_scale = st.slider('BondScale',min_value=0.0,max_value=1.0,value=0.8)
        brightness = st.slider('Brightness',min_value=0.0,max_value=1.0,value=0.4)
        relativeAtomScale = st.slider('Relative atom scale',min_value=0.0,max_value=1.0,value=0.64)
        atomShade  = st.slider('Atom shade',min_value=0.0,max_value=1.0,value=0.5)

    _PARAMETERS = {'outline': outl , 'bondScale': bond_scale,
                    'bonds': bond ,'atomShade' : atomShade,
                    'brightness': brightness, 'relativeAtomScale':relativeAtomScale,
                    }
    res = spec_plot(example_xyz,_PARAMETERS = _PARAMETERS,wbox_height="500px",wbox_width="500px")
        '''
    )
