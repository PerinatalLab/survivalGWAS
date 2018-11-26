#!/usr/bin/python3

# This script is used to format summary files obtained from the script cox_parallel.R
# SNP id will be defined as:
# CHR:POS:EFF:REF with EFF alphabetically higher than REF

import pandas as pd

def digest ( df ):


