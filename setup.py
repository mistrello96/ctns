from setuptools import setup, find_packages
import os

with open('requirements.txt') as f:
    required = f.read().splitlines()

with open('README.md') as readme_file:
    README = readme_file.read()

ver = os.environ['CI_COMMIT_TAG']

setup(
  name = 'ctns',
  packages = find_packages(),
  scripts = ['bin/ctns', 'bin/ctns.bat','bin/ctns_run.py'],
  description = 'CTNS, Contact Tracing Network Simulator: a tool to simulate contact tracing in a population where a disease is spreading',
  long_description_content_type = "text/markdown",
  long_description = README,
  license = 'GNU GPLv3',
  author = 'Matteo Mistri, Diego Miglio',
  author_email = 'matteo.mistri1996@gmail.com',
  install_requires = required,
  version = ver,
  url = "https://gitlab.com/mistrello96/ctns",
  download_url='https://pypi.org/project/ctns/',
)