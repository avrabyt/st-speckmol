import streamlit as st
from PIL import Image
st.markdown('''# st-speckmol :package:
_A Streamlit **Component** for creating Speck molecular structures within Streamlit Web app._
''')
image = Image.open("Resources/mol.png")
st.image(image)
st.sidebar.header(":sparkles: Install")
st.sidebar.code('''pip install st-speckmol''')
st.sidebar.subheader(':rocket: Usage')
st.sidebar.code('from st_speckmol import spec_plot')
st.sidebar.markdown(":heavy_check_mark: Quick Note: The module is **'spec_plot'** not ~~spec**k**_plot~~")
st.sidebar.markdown(
    "[![](https://static.pepy.tech/badge/st-speckmol)](https://pypi.org/project/st-speckmol/) "

    " [![](https://static.pepy.tech/badge/st-speckmol/month)](https://pepy.tech/project/st-speckmol)"
    
)
st.sidebar.markdown(
    '''
    [![Follow](https://img.shields.io/twitter/follow/Avra_b?style=social)](https://www.twitter.com/Avra_b)
    [![Buy me a coffee](https://img.shields.io/badge/Buy%20me%20a%20coffee--yellow.svg?logo=buy-me-a-coffee&logoColor=orange&style=social)](https://www.buymeacoffee.com/AvraCodes) 
    '''
)

st.info(":bug: :hankey:  Report any bugs/issues here - [st-speckmol](https://github.com/avrabyt/st-speckmol)")
