# set_idd.py
import eppy
from eppy import modeleditor
from eppy.modeleditor import IDF

# Define the IDD file path using raw strings
iddfile = r'C:\Users\alienor\Documents\hess\resources\loads\Energyplus.idd'

# Set the IDD file
IDF.setiddname(iddfile)

print(f"IDD file has been set to: {iddfile}")
