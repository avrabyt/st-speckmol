import setuptools

# read the contents of your README file
with open("README.md", "r", encoding='utf8') as fh:
    long_description = fh.read()

setuptools.setup(
    name="st_speckmol",
    version="0.0.5",
    author="Avratanu Biswas",
    author_email="avrab.yt@gmail.com",
    description="Streamlit component for for Speck molecule visualization.",
    long_description=long_description,
    long_description_content_type='text/markdown',
    url="https://github.com/avrabyt/Specklit",
    packages=setuptools.find_packages(),
    include_package_data=True,
    classifiers=[],
    python_requires=">=3.6",
    install_requires=["streamlit >= 0.63", "ipyspeck==0.6.1", "ipywidgets==7.6.3"],
)
