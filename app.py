import streamlit as st
import ipywidgets as widgets
from ipywidgets import embed
import ipyspeck
import streamlit.components.v1 as components
import glob

st.set_page_config(
    layout="centered",
    page_title="Specklit",
    page_icon=":sparkles:")

def add_spec_param(spec_xyz,atom_Scale,bond_Scale,outline,atomShade,bonds,bondThreshold):
    #Modify atoms size
    spec_xyz.atomScale = atom_Scale
    #change bonds size
    spec_xyz.bondScale = bond_Scale
    #highlight borders
    spec_xyz.outline = outline
    spec_xyz.atomShade = atomShade
    spec_xyz.bonds = bonds
    spec_xyz.bondThreshold = bondThreshold
    return spec_xyz

def spec_plot(_xyz):
    spec_xyz = ipyspeck.speck.Speck(data=_xyz, title='')
    spec_xyz = add_spec_param(spec_xyz,atom_Scale,bond_Scale,outline,atom_Shade,bonds,bondThreshold)
    spec_xyz.setColorSchema(schema='newcpk')
    c = widgets.Box([spec_xyz], layout=widgets.Layout(width="800px",height="700px"))
    snippet = embed.embed_snippet(c)
    html = embed.html_template.format(title="", snippet=snippet)
    components.html(html,height = 1000, width = 900)
    return spec_xyz

ex_files = glob.glob("examples/*.xyz")

st.sidebar.title(" XYZ to SPECK-Structures :tada:")
with st.sidebar.expander(label = "Examples",expanded=True):
    st.markdown("[Source](https://github.com/wwwtyro/speck/tree/gh-pages/static/samples)")
    example_xyz = st.selectbox("Select a molecule",ex_files)
    f = open(example_xyz,"r")
    example_xyz = f.read()


st.sidebar.markdown("# Parameter section:")
st.sidebar.info("You can also add your own coordinates below. :coffee:")
_xyz = st.sidebar.text_area(
                label = "What are the Coordinates ?",
                value= example_xyz, height  = 200)

st.code(_xyz.splitlines()[1])


atom_Scale = st.sidebar.slider('Atom Scale', 0.1, 1.0, 0.35)
bond_Scale = st.sidebar.slider('Bond Scale', 0.1, 1.0, 0.5)
atom_Shade = st.sidebar.slider('Atom Shade', 0.1, 1.0, 0.0)
bondThreshold = st.sidebar.slider('Bond Threshold', 0.1, 5.0, 1.2)
bonds = st.sidebar.checkbox("Bonds",value = True)
outline = st.sidebar.checkbox('Outline',value = True)
res = spec_plot(_xyz)

with st.expander("References:",expanded=True):
    st.markdown("[Speck Online](http://wwwtyro.github.io/speck/)")
    st.markdown("[Speck Python package](https://pypi.org/project/ipyspeck/)")