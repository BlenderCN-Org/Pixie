# Pixie
Experimental renderer for FX.

This project is a learning exercise as much as anything else so everything has been done from scratch where possible (with some exceptions such as numpy to export binary files and pillow to save images to disk).

NOTE: All files paths are hard coded in the python files so you will need to set those before starting.

To use load smokeExporter.py into Blender with a scene containing a domain object, select that domain and run the script. Then (this has only been tested on linux) run this command in a terminal:
python3 setup.py build && python3 testModule.py

