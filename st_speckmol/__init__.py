import streamlit.components.v1 as components
import ipywidgets as widgets
from ipywidgets import embed
import ipyspeck

# # DEFAULT VALUES 
# atom_Scale = 0.5
# bond_Scale = 0.5
# outline = False
# atomShade = 0.0
# bonds = True
# bondThreshold = 1.2

# def add_speck_param(spec_xyz,atom_Scale,bond_Scale,outline,atomShade,bonds,bondThreshold):
    
#     #Modify atoms size
#     spec_xyz.atomScale = atom_Scale
#     #change bonds size
#     spec_xyz.bondScale = bond_Scale
#     #highlight borders
#     spec_xyz.outline = outline
#     spec_xyz.atomShade = atomShade
#     spec_xyz.bonds = bonds
#     spec_xyz.bondThreshold = bondThreshold
#     return spec_xyz

class st_speck:
    # # # DEFAULT VALUES 
    # atom_Scale = 0.5
    # bond_Scale = 0.5
    # outline = False
    # atomShade = 0.0
    # bonds = True
    # bondThreshold = 1.2
    # specobj = ipyspeck.speck.Speck()
    def speck_plot(self,_xyz, wbox_height="700px", 
                wbox_width="800px",
                component_h = 700, 
                component_w = 800, 
                scroll = False):

        """ Plots the speckmol molecule using the ipyspeck library and returns 
            the <class 'ipyspeck.speck.Speck'>.
        
        Parameters
        ----------
        _xyz : str
            The xyz string of the molecule.
        wbox_height : str
            The height of the widget box.
        wbox_width : str
            The width of the widget box.
        component_h : int
            The height of the streamlit html component.
        component_w : int
            The width of the streamlit html component.
        scroll : bool
            If True, the streamlit component will scroll.  

        Returns
        -------
        spec_xyz : <class 'ipyspeck.speck.Speck'>
            The speckmol molecule. spec_xyz.keys() returns the keys of the
            molecule. For example - spec_xyz.keys() returns ['atomScale',  
            'bondScale', 'atomShade', 'bondThreshold', 'bondColor', 'atomColor',
            'outline', 'bonds', 'atomScale', 'atomColor', 'atomScale', 'atomScale]
            These keys are useful for modifying the molecule.

        """        
        # Read the xyz file
        spec_xyz = ipyspeck.speck.Speck(data = _xyz)
        # Add Parameters
        # add_spec_param(spec_xyz,atom_Scale,bond_Scale,outline,atomShade,bonds,bondThreshold)
        # Create the widget box
        widg = widgets.Box([spec_xyz], layout=widgets.Layout(height=wbox_height,width=wbox_width))
        # Embed the widget box in the streamlit html component
        sc = embed.embed_snippet(widg)
        html = embed.html_template.format(title="", snippet=sc)
        components.html(html,height = component_h, width = component_w,scrolling=scroll)
        return spec_xyz

    def add_speck_param(self):
        speckobj = ipyspeck.speck.Speck()
        #Modify atoms size
        # specobj.atomScale = atom_Scale
        # #change bonds size
        # specobj.bondScale = bond_Scale
        # #highlight borders
        # specobj.outline = outline
        # specobj.atomShade = atomShade
        # specobj.bonds = bonds
        speckobj.bondThreshold = 1.6
            # return self


import glob
import streamlit as st
ex_files = glob.glob("examples/*.xyz")
with st.sidebar:
    example_xyz = st.selectbox("Select a molecule",ex_files)
    f = open(example_xyz,"r")
    example_xyz = f.read()

# res = speck.speck_plot(example_xyz)
# st.write(res.keys)

# z = st_speck()
# st.write(z.atomScale)
# z.add_speck_param(z)
# z.bondThreshold

# res = z.speck_plot(example_xyz)
# res.bondThreshold

z = st_speck()
z.add_speck_param()
st.write(z.add_speck_param)

# st.write(z.speck_plot(example_xyz))