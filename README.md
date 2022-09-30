# Stspeckmol
A Streamlit Component for creating Speck molecular structures within Streamlit Web app.

<table>
    <tr>
        <td>Latest Release</td>
        <td>
            <a href="https://pypi.org/project/st-speckmol/"/>
            <img src="https://static.pepy.tech/badge/st-speckmol"/>
        </td>
    </tr>
    <tr>
        <td>PyPI Downloads</td>
        <td>
            <a href="https://pepy.tech/project/st-speckmol"/>
            <img src="https://static.pepy.tech/badge/st-speckmol/month"/>
        </td>
    </tr>
</table>

## Installation 
```console
pip install st-speckmol
```
to upgrade use,
```console
pip install --upgrade st-speckmol
```

## Example

Try the app, for different examples. 

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://share.streamlit.io/avrabyt/specklit/main/app.py)


## Quickstart

``` python
import streamlit as st
import glob
from st_speckmol import speck_plot

# Example files path
ex_files = glob.glob("examples/*.xyz")
with st.sidebar:
    example_xyz = st.selectbox("Select a molecule",ex_files)
    f = open(example_xyz,"r")
    example_xyz = f.read()

res = speck_plot(example_xyz)

```

![Speclit demo](https://github.com/avrabyt/Specklit/blob/main/Resources/SpeckLit_demo.gif)

## Video tutorial
<details>
  <summary>See the tutorial video</summary>

[![How to Build PROTEIN VISUALIZATION WEB-APP using PYTHON and STREAMLIT | PART 1](https://github.com/avrabyt/st-speckmol/blob/main/Resources/Speck-Thumbnail.png)](https://youtu.be/jUh923Z4fuk)

</details>


------------------------

## Related library - [Stmol](https://github.com/napoles-uach/stmol) 

During the development of the related and popular library [Stmol](https://github.com/napoles-uach/stmol), we introduced `speck_plot()` function for easy usage of both libraries simultaneously. However, the entire StSpeckmol has not yet been merged and extra(read new) functions such as `add_speck_param` is only available with this library.

```python
# Installation of Stmol
pip install stmol==0.0.9

# Import Speck plot
from stmol import speck_plot

```

Incase you are using `StSpeckmol` for scientific purposes for speck visualization, make sure you use `Stmol` (https://doi.org/10.3389/fmolb.2022.990846) and cite as following, 
```console
Nápoles-Duarte JM, Biswas A,Parker MI, Palomares-Baez JP, Chávez-Rojo MA and Rodríguez-Valdez LM (2022), 
Stmol: A component for building interactive molecular visualizations within streamlit web-applications.
Front. Mol. Biosci. 9:990846. doi: 10.3389/fmolb.2022.990846
```





_Note : Meanwhile,[ipyspeck](https://pypi.org/project/ipyspeck/) in their latest release ` 0.6.1 ` has added the stspec module [https://github.com/avrabyt/Specklit/issues/1#issuecomment-1134798584], therefore feel free to use whatever convinient, as long as you are interested to have fun with beautiful speck strcutures 🧬 and streamlit_ 🎈 🎉 


------------------------
## References

[Speck Online](http://wwwtyro.github.io/speck/)

[Speck Python package](https://pypi.org/project/ipyspeck/)

[Example-Source](https://github.com/wwwtyro/speck/tree/gh-pages/static/samples)

[Stmol](https://github.com/napoles-uach/stmol)

