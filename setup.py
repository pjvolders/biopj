from setuptools import setup, find_packages

setup(name='biopj',
      version='0.1',
      description='The funniest joke in the world',
      # url='http://github.com/storborg/funniest',
      author='Pieter-Jan Volders',
      author_email='pieterjan.volders@ugent.be',
      license='MIT',
      packages=find_packages(),
      entry_points={
            'console_scripts': [
                  'circ_seq = biopj.circ_seq:main'
            ]
      },
      zip_safe=False)