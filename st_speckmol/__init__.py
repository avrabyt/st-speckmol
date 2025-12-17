from .utils import *
import streamlit.components.v1 as components
from .version import __version__

# Try to import the native Streamlit component from ipyspeck
try:
    from ipyspeck.stspeck import Speck as IpySpeck
    _HAS_STSPECK = True
except ImportError:
    _HAS_STSPECK = False
    # Fallback to ipywidgets approach (works locally but not on Streamlit Cloud)
    import ipywidgets as widgets
    from ipywidgets import embed
    import ipyspeck

def speck_plot(_xyz, wbox_height="700px", 
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
    
    res = speck_plot(example_xyz,_PARAMETERS = _PARAMETERS,wbox_height="500px",wbox_width="500px") # with _PARAMETERS
    
    res = speck_plot(example_xyz) # without _PARAMETERS or Default Parameters
    
    Returns
    -------
    speck_xyz : <class 'ipyspeck.speck.Speck'>
        The speckmol molecule. speck_xyz.keys() returns the keys of the
        molecule. For example - speck_xyz.keys() returns ['atomScale',  
        'bondScale', 'atomShade', 'bondThreshold', 'bondColor', 'atomColor',
        'outline', 'bonds', 'atomScale', 'atomColor', 'atomScale', 'atomScale]
        These keys are useful for modifying the molecule.
    
    """        
    # Use native Streamlit component if available (works on Streamlit Cloud)
    if _HAS_STSPECK:
        # Extract parameters from _PARAMETERS dict if provided
        params = {}
        if _PARAMETERS:
            params = dict(*_PARAMETERS.values())
        
        # Set dimensions
        params['width'] = wbox_width
        params['height'] = wbox_height
        
        # Call the native Streamlit component
        IpySpeck(data=_xyz, **params)
        
        # Return a mock object for compatibility (native component doesn't return speck object)
        class MockSpeck:
            def __init__(self):
                self.data = _xyz
        return MockSpeck()
    else:
        # Fallback to ipywidgets approach (works locally but not on Streamlit Cloud)
        speck_xyz = ipyspeck.speck.Speck(data = _xyz) 
        if _PARAMETERS :
            param = dict(*_PARAMETERS.values())
            add_speck_param(speck_xyz,kwargs = param)
        # Create the widget box
        widg = widgets.Box([speck_xyz], layout=widgets.Layout(height=wbox_height,width=wbox_width))
        # Embed the widget box in the streamlit html component
        sc = embed.embed_snippet(widg)
        html = embed.html_template.format(title="", snippet=sc)
        components.html(html,height = component_h, width = component_w,scrolling=scroll)
        return speck_xyz

def add_speck_param(xyz, **kwargs):
    '''
    Helper function for speck_plot. Called within the function to add
    parameters / attributes to the molecules.
    
    Attributes
    ----------
    _view_name : str
        (DOMWidget) Name of the widget view class in front-end
    _model_name :dofStrength
        Name of the widget model class in front-end
    _view_module : str
        (DOMWidget) Name of the front-end module containing widget view
    _model_module : str
        (DOMWidget) Name of the front-end module containing widget model
    _view_module_version : str
        (DOMWidget) Version of the front-end module containing widget view
    _model_module_version : str
        (DOMWidget) Version of the front-end module containing widget model
    data : str
        xyz model, default(True)
    bonds : bool
        Enable visualizations of bonds?, default(True)
    atomScale : float
        Atom radius, size of spheres, default(0.24)
    relativeAtomScale : float
        Relative atom radius, default(0.64)
    bondScale : float
        bonds size, size of the tubes connecting atoms, default(0.5)
    brightness : float
        brightness, default(0.5)
    outline : float
        Outline strength, default(0.0)
    spf : float
        Samples per frame, default(32)
    bondThreshold : float
        Bonding radius, defines the max distance for atoms to be connected,
        default(1.2)
    bondShade : float
        bonds shade, default(0.5)
    atomShade : float
        Atoms shade, default(0.5)
    dofStrength : float
        Depth of field strength, default(0.0)
    dofPosition : float
        Depth of field position, default(0.5)

    Reference - https://github.com/denphi/speck/blob/master/widget/ipyspeck/ipyspeck/speck.py    

    '''
    param = dict(*kwargs.values())
    for key,value in param.items():
        if key == 'outline':
            xyz.outline = value 
        if key == 'bonds':
            xyz.bonds = value   
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



@deprecated('You must use speck_plot instead. spec_plot will be deprecated in future release. ')
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
    # Use native Streamlit component if available (works on Streamlit Cloud)
    if _HAS_STSPECK:
        # Extract parameters from _PARAMETERS dict if provided
        params = {}
        if _PARAMETERS:
            params = dict(*_PARAMETERS.values())
        
        # Set dimensions
        params['width'] = wbox_width
        params['height'] = wbox_height
        
        # Call the native Streamlit component
        IpySpeck(data=_xyz, **params)
        
        # Return a mock object for compatibility
        class MockSpeck:
            def __init__(self):
                self.data = _xyz
        return MockSpeck()
    else:
        # Fallback to ipywidgets approach (works locally but not on Streamlit Cloud)
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

@deprecated('You must use add_speck_param instead. add_spec_param will be deprecated in future release. ')
def add_spec_param(xyz, **kwargs):
    '''
    Helper function for spec_plot. Called within the function to add
    parameters / attributes to the molecules.
    
    Attributes
    ----------
    _view_name : str
        (DOMWidget) Name of the widget view class in front-end
    _model_name :dofStrength
        Name of the widget model class in front-end
    _view_module : str
        (DOMWidget) Name of the front-end module containing widget view
    _model_module : str
        (DOMWidget) Name of the front-end module containing widget model
    _view_module_version : str
        (DOMWidget) Version of the front-end module containing widget view
    _model_module_version : str
        (DOMWidget) Version of the front-end module containing widget model
    data : str
        xyz model, default(True)
    bonds : bool
        Enable visualizations of bonds?, default(True)
    atomScale : float
        Atom radius, size of spheres, default(0.24)
    relativeAtomScale : float
        Relative atom radius, default(0.64)
    bondScale : float
        bonds size, size of the tubes connecting atoms, default(0.5)
    brightness : float
        brightness, default(0.5)
    outline : float
        Outline strength, default(0.0)
    spf : float
        Samples per frame, default(32)
    bondThreshold : float
        Bonding radius, defines the max distance for atoms to be connected,
        default(1.2)
    bondShade : float
        bonds shade, default(0.5)
    atomShade : float
        Atoms shade, default(0.5)
    dofStrength : float
        Depth of field strength, default(0.0)
    dofPosition : float
        Depth of field position, default(0.5)

    Reference - https://github.com/denphi/speck/blob/master/widget/ipyspeck/ipyspeck/speck.py    

    '''
    param = dict(*kwargs.values())
    for key,value in param.items():
        if key == 'outline':
            xyz.outline = value 
        if key == 'bonds':
            xyz.bonds = value   
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

