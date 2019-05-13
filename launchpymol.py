# -*- coding: utf-8 -*-
"""
Created on Wed Sep 26 00:44:56 2012

@author: Afsar
"""

#!/usr/bin/env python
 
# Tell PyMOL we don't want any GUI features.
import __main__
__main__.pymol_argv = [ 'pymol','-r', 'pymoltest.py']
 
# Importing the PyMOL module will create the window.
import pymol
 
# Call the function below before using any PyMOL modules.
pymol.finish_launching()
 
from pymol import cmd

cmd.load("Ltemp.pdb")
cmd.load("Rtemp.pdb")

cmd.show('cartoon',"all")

