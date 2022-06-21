import streamlit.components.v1 as components
import ipywidgets as widgets
from ipywidgets import embed
import ipyspeck

def spec_plot(_xyz, wbox_height="700px", 
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
    # Create the widget box
    widg = widgets.Box([spec_xyz], layout=widgets.Layout(height=wbox_height,width=wbox_width))
    # Embed the widget box in the streamlit html component
    sc = embed.embed_snippet(widg)
    html = embed.html_template.format(title="", snippet=sc)
    components.html(html,height = component_h, width = component_w,scrolling=scroll)
    return spec_xyz