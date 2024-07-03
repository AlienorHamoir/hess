# set_idd.py
import eppy
from eppy import modeleditor
from eppy.modeleditor import IDF

# Define the IDD file path using raw strings
iddfile = r'/home/alienor/Documents/hess/resources/loads/V1-4-0-Energy+.idd'

# Set the IDD file
IDF.setiddname(iddfile)

print(f"IDD file has been set to: {iddfile}")
