from setuptools import setup, find_packages

setup(
    name="cablab-oz",
    version="0.1.0",
    description='CAB-LAB Ozone Extension',
    license='THE BEER-WARE LICENSE',
    author='John Doe',
    packages=find_packages(),
    entry_points={
        'cablab.converter.classes': [
            'oz = my_cablab_ext:OzoneConverter',
        ],
    }
)
