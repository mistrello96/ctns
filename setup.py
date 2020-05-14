from setuptools import setup, find_packages
import os

with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
  name = 'ctns',
  packages = find_packages(),
  version = '0.2',
  scripts = ['bin/ctns', 'bin/ctns.bat','bin/ctns_run.py'],
  description = 'CTNS, Contact Tracing Network Simulator: a tool to simulate contact tracing in a population where a disease is spreading',
  author = 'Matteo Mistri, Diego Miglio',
  author_email = 'matteo.mistri1996@gmail.com',
  install_requires = required,
  url = "https://gitlab.com/mistrello96/ctns",
  download_url='https://pypi.org/project/ctns/',
)