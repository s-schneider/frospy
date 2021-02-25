rm -rf dist build
python setup.py sdist bdist_wheel

twine check dist/*
twine upload --repository-url https://test.pypi.org/legacy/ dist/*



# for django
conda create --name django python=3.8
conda activate django
conda config --add channels conda-forge
conda install pip
conda install obspy
conda install cartopy
conda install basemap


python setup.py install

pip install psycopg2-binary
pip install pyshtools
