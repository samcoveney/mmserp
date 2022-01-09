#from setuptools import setup, find_packages
from distutils.core import setup

setup(name = 'mmserp',
      version = '1.0',
      description = 'Tools for ERP inference with the modified Mitchell-Schaeffer model',
      author = 'Sam Coveney',
      author_email = 'coveney.sam@gmail.com',
      license = 'GPL-3.0+',
      packages = ['mmserp'],
      package_dir = {'mmserp': 'mmserp'},
      package_data = {'mmserp': ['data/surrogate_coefficients.npz', 'data/stanmodel_tophat.pkl', 'data/stanmodel_gaussian.pkl', 'data/runCARP.sh']},
      scripts=['scripts/mmserp_viewMesh', 'scripts/mmserp_viewEigs', 'scripts/mmserp_createCARPfiles', 'scripts/mmserp_generateFields', 'scripts/mmserp_browseHDF5', 'scripts/mmserp_meshToHDF5', 'scripts/mmserp_duplicateHDF5', 'scripts/mmserp_createStimulus', 'scripts/mmserp_decimateMesh', 'scripts/mmserp_getSimResults', 'scripts/mmserp_inference'],
     )

