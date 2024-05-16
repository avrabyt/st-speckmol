# Stspeckmol
A Streamlit Component for creating Speck molecular structures within Streamlit Web app.

[![PyPI version](https://badge.fury.io/py/st-speckmol.svg)](https://pypi.org/project/st-speckmol/)
[![Downloads](https://pepy.tech/badge/st-speckmol)](https://pepy.tech/project/st-speckmol)
[![Downloads](https://pepy.tech/badge/st-speckmol/month)](https://pepy.tech/project/st-speckmol)
[![Website shields.io](https://img.shields.io/website-up-down-green-red/http/shields.io.svg)](https://hellostspeckmol.streamlitapp.com)
![GitHub last commit](https://img.shields.io/github/last-commit/avrabyt/st-speckmol?style=plastic)
![GitHub Release Date](https://img.shields.io/github/release-date/avrabyt/st-speckmol)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11205344.svg)](https://doi.org/10.5281/zenodo.11205344)



## Installation 
```console
pip install st-speckmol
```
to upgrade use,
```console
pip install --upgrade st-speckmol
```
>:warning: https://github.com/avrabyt/st-speckmol/issues/20 In case of `ModuleNotFoundError: No module named 'ipython_genutils'` :
>```console
> pip install ipython_genutils
>```
> **Future release of st-speckmol, will support this module natively.**

## Example

Try the app, for different examples. 

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://hellostspeckmol.streamlitapp.com)


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

[Blog Post](https://medium.com/@avra42/how-to-build-molecular-structures-visualizing-web-application-using-python-and-streamlit-5ec9da86550c) 


## Scientfic usage
To cite any scientific usage, please refer to the following citation:

```
@software{Biswas_2024,
  author = {Avratanu Biswas},
  title = {st-speckmol},
  version = {v0.0.6.1},
  doi = {10.5281/zenodo.11205344},
  url = {https://github.com/avrabyt/st-speckmol},
  date = {2024-05-16}
}

```
For additional assistance, feel free to reach out to me directly.

## Related library - [Stmol](https://github.com/napoles-uach/stmol) 
During the development of the related library [Stmol](https://github.com/napoles-uach/stmol), we introduced the convenient `speck_plot()` function, allowing seamless integration of both libraries. 
> ℹ️
> The complete integration of StSpeckmol **has not yet been finalized**, and any additional (or recently implemented) functionalities like `add_speck_param` will exclusively be accessible through this library."

```python
# Installation of Stmol
pip install stmol==0.0.9

# Import Speck plot
from stmol import speck_plot

```

## References

[Speck Online](http://wwwtyro.github.io/speck/)

[Speck Python package](https://pypi.org/project/ipyspeck/)

[Example-Source](https://github.com/wwwtyro/speck/tree/gh-pages/static/samples)

[Stmol](https://github.com/napoles-uach/stmol)

