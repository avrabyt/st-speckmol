import py3Dmol
import streamlit as st


def visbox2(objeto,bxi,byi,bzi,bxf,byf,bzf):  
    objeto.addBox({'center':{'x':bxi,'y':byi,'z':bzi},'dimensions': {'w':bxf,'h':byf,'d':bzf},'color':'blue','opacity': 0.5})

def complxvis(objeto,protein_name,ligand_name):
    mol1 = open(protein_name, 'r').read()
    mol2 = open(ligand_name, 'r').read()
    objeto.addModel(mol1,'pdb')
    objeto.setStyle({'cartoon': {'color':'spectrum'}})
    objeto.addModel(mol2,'pdb')
    objeto.setStyle({'model':1},{'stick':{}})

def vismol(bxi=-10,byi=-10,bzi=-10,bxf=5,byf=5,bzf=5):  
    mol_view = py3Dmol.view(800, 400,viewergrid=(1,2))  
    visbox2(mol_view,bxi,byi,bzi,bxf,byf,bzf)
    complxvis(mol_view,'3lfa.pdb','ligand.pdb')
    mol_view.setBackgroundColor('0xeeeeee')
    mol_view.rotate(90, {'x':0,'y':1,'z':0},viewer=(0,1));
    mol_view.zoomTo()  
    mol_view.show()

from ipywidgets import interact,fixed,IntSlider,embed,widgets
import ipywidgets
import streamlit.components.v1 as components
widg = interact(vismol,bxi=ipywidgets.IntSlider(min=-100,max=100, step=1) ,byi=ipywidgets.IntSlider(min=-100,max=100, step=1),
        bzi=ipywidgets.IntSlider(min=-100,max=100, step=1),bxf=ipywidgets.IntSlider(min=0,max=30, step=1),
        byf=ipywidgets.IntSlider(min=0,max=30, step=1),
        bzf=ipywidgets.IntSlider(min=0,max=30, step=1))
st.write(widg)
# widg = widgets.Box([cc], layout=widgets.Layout(height="700px",width="700px"))
sc = embed.embed_snippet(widg)
html = embed.html_template.format(title="", snippet=sc)
components.html(html)