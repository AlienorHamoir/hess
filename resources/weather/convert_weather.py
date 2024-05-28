import os
import sys
import subprocess as sp

root = os.path.dirname(os.path.abspath(__file__))
for f in os.listdir():
    if f.endswith('epw'):
        sp.call('java -jar ../modelica/Buildings/Resources/bin/ConvertWeatherData.jar {}'.format(f), shell=True)

