# Specklit
Speck figures to Streamlit Web App

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://share.streamlit.io/avrabyt/specklit/main/app.py)

![Speclit demo](https://github.com/avrabyt/Specklit/blob/main/SpeckLit_demo.gif)

## Installation
`pip install st-speckmol==0.0.3`

## Example

```
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


# References

[Speck Online](http://wwwtyro.github.io/speck/)

[Speck Python package](https://pypi.org/project/ipyspeck/)

[Example-Source](https://github.com/wwwtyro/speck/tree/gh-pages/static/samples)
