# Specklit
Speck figures to Streamlit Web App

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
pip install st-speckmol==0.0.3
```

## Example

``` python
import streamlit as st
import glob
from st_speckmol import spec_plot

# Example files path
ex_files = glob.glob("examples/*.xyz")
with st.sidebar:
    example_xyz = st.selectbox("Select a molecule",ex_files)
    f = open(example_xyz,"r")
    example_xyz = f.read()

res = spec_plot(example_xyz)

```
[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://share.streamlit.io/avrabyt/specklit/main/app.py)

![Speclit demo](https://github.com/avrabyt/Specklit/blob/main/Resources/SpeckLit_demo.gif)

Meanwhile,[ipyspeck](https://pypi.org/project/ipyspeck/) in their latest release ` 0.6.1 ` has added the stspec module [https://github.com/avrabyt/Specklit/issues/1#issuecomment-1134798584], therefore feel free to use whatever convinient, as long as you are interested to have fun with beautiful speck strcutures 🧬 and streamlit 🎈 🎉

# References

[Speck Online](http://wwwtyro.github.io/speck/)

[Speck Python package](https://pypi.org/project/ipyspeck/)

[Example-Source](https://github.com/wwwtyro/speck/tree/gh-pages/static/samples)

