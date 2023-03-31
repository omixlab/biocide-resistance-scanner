from setuptools import setup, find_packages

setup(
    name="brs",
    version='0.0.15',
    packages=find_packages(),
    author="Elias Eduardo Barbosa, Frederico Schmitt Kremer",
    author_email="fred.s.kremer@gmail.com",
    description="Biocide Resistance Scanner",
    long_description=open("README.md").read(),
    long_description_content_type='text/markdown',
    keywords="bioinformatics machine-learning data science drug discovery QSAR",
    entry_points = {'console_scripts':[
        'brs = brs.main:main',
        ]},
    install_requires = [
        requirement.strip('\n') for requirement in open("requirements.txt")
    ]
)