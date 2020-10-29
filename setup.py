from setuptools import setup, find_packages

setup(name='biopj',
      version='0.1',
      description='Some useful homegrown Python modules and scripts',
      author='Pieter-Jan Volders',
      author_email='pieterjan.volders@ugent.be',
      license='MIT',
      packages=find_packages(),
      entry_points={
            'console_scripts': [
                  'circ_seq = biopj.circ_seq:main',
                  'reverse_complement = biopj.reverse_complement:main',
                  'translate = biopj.translate:main',
                  'split_at_stop_codons = biopj.proteomics.split_at_stop_codons:main'
            ]
      },
      zip_safe=False)