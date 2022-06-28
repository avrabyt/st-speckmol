import requests
#import openbabel
import streamlit.components.v1 as components
import ipywidgets as widgets
from ipywidgets import embed
import ipyspeck


def pdb2xyz(pdb_file,type='.pdb'):
    """
    Reads a pdb file and returns the xyz string.
    """
    if type == '.pdb':
        url = "https://files.rcsb.org/view/" + pdb_file + type
    elif type == '.sdf':
        url = "https://files.rcsb.org/ligands/view/" + pdb_file + type 
    else:
        None
    
    if not(None):       
        r = requests.get(url)
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("pdb", "xyz")
        mol = openbabel.OBMol()
        obConversion.ReadString(mol, r.text)
        xyz_mol = obConversion.WriteString(mol)
    return  xyz_mol

def add_spec_param(xyz, **kwargs):
    param = dict(*kwargs.values())
    for key,value in param.items():
        if key == 'outline':
            xyz.outline = value 
        if key == 'bonds':
            xyz.bonds = value  # Bool (0|1) 
        if key == 'atomScale':
            xyz.atomScale = value
        if key == 'relativeAtomScale':
            xyz.relativeAtomScale = value
        if key == 'bondScale':
            xyz.bondScale = value
        if key == 'brightness':
            xyz.brightness = value
        if key == 'spf':
            xyz.spf = value
        if key == 'bondThreshold':
            xyz.bondThreshold = value
        if key == 'bondShade':
            xyz.bondShade = value
        if key == 'atomShade':
            xyz.atomShade = value
        if key == 'dofStrength':
            xyz.dofStrength = value
        if key == 'dofPosition':
            xyz.dofPosition = value    
        if key == 'width':
            xyz.width = value
        if key == 'height':
            xyz.height = value
    return None

def spec_plot(_xyz, wbox_height="700px", 
            wbox_width="800px",
            component_h = 700, 
            component_w = 800, 
            scroll = False, 
            **_PARAMETERS : dict):

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
    _PARAMETERS : Dict
        Define custom Parameters within a dictionary named "_PARAMETERS"
        
    Example
    ---------
    _PARAMETERS = {'outline': True , 'bondScale': 0.9, 'bonds': False}
    
    res = spec_plot(example_xyz,_PARAMETERS = _PARAMETERS,wbox_height="500px",wbox_width="500px") # with _PARAMETERS
    
    res = spec_plot(example_xyz) # without _PARAMETERS or Default Parameters
    
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
    if _PARAMETERS :
        param = dict(*_PARAMETERS.values())
        add_spec_param(spec_xyz,kwargs = param)
    # Create the widget box
    widg = widgets.Box([spec_xyz], layout=widgets.Layout(height=wbox_height,width=wbox_width))
    # Embed the widget box in the streamlit html component
    sc = embed.embed_snippet(widg)
    html = embed.html_template.format(title="", snippet=sc)
    components.html(html,height = component_h, width = component_w,scrolling=scroll)
    return spec_xyz