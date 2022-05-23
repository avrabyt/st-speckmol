import setuptools

setuptools.setup(
    name="st_speckmol",
    version="0.0.3",
    author="Avratanu Biswas",
    author_email="avrab.yt@gmail.com",
    description="Streamlit component for for Speck molecule visualization.",
    long_description="Streamlit component for visualizing SPECK molecules - https://github.com/wwwtyro/speck",
    long_description_content_type="text/plain",
    url="https://github.com/avrabyt/Specklit",
    packages=setuptools.find_packages(),
    include_package_data=True,
    classifiers=[],
    python_requires=">=3.6",
    install_requires=["streamlit >= 0.63", "ipyspeck==0.6.1", "ipywidgets==7.6.3"],
)
