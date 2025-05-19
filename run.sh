#!/bin/zsh

# Get all the data for the Be7 CN
python3 fetch_exfor.py Be7 --EnergyMin 0.01 --EnergyMax 15 --pops ./ripl2pops.xml --dir data/X4

# Convert to AZURE2 format
python3 prepare_azure.py

# Convert to SAMMY format
python3 prepare_sammy.py