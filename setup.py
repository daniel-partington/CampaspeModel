try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup
    import re
    import os

    def find_packages(path='.'):
        ret = []
        for root, dirs, files in os.walk(path):
            if '__init__.py' in files:
                ret.append(re.sub('^[^A-z0-9_]+', '', root.replace('/', '.')))
        return ret
    # End find_packages()
# End try

config = {
    'description': 'CampaspeModel',
    'author': 'Daniel Partington, Takuya Iwanaga',
    'url': 'https://github.com/daniel-partington/CampaspeModel',
    'download_url': 'https://github.com/daniel-partington/CampaspeModel/archive/master.zip',
    'author_email': 'dpartington1982@gmail.com',
    'version': '0.1',
    'install_requires': ['flopy', 'numpy', 'pandas', 'gdal', 'matplotlib'],
    'packages': find_packages(),
    'scripts': [],
    'name': 'CampaspeModel'
}

setup(**config)
