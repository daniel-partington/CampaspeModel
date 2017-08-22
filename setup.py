try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'CampaspeModel',
    'author': 'Daniel Partington, Takuya Iwanaga',
    'url': 'https://github.com/daniel-partington/CampaspeModel',
    'download_url': 'https://github.com/daniel-partington/CampaspeModel/archive/master.zip',
    'author_email': 'dpartington1982@gmail.com',
    'version': '0.1',
    'install_requires':['flopy', 'numpy', 'pandas', 'gdal', 'matplotlib'],
    'packages': ['CampaspeModel', 'CampaspeModel.GW_link_Integrated'],
    'scripts': [],
    'name': 'CampaspeModel'
    }
	
setup(**config)	