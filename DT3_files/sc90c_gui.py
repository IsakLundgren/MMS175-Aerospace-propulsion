import os
from pathlib import Path

sc90c_working_dir = Path.cwd() / 'sc90c_files'
os.chdir(sc90c_working_dir)
## Launches the extra command window that is needed for createInputSC90C.py
os.system('start "Command Prompt" cmd.exe 1')
## Launches the python script createInputSC90C.py
os.system('start "Controls createInputSC90C.py" python createInputSC90C.py 2')
##Launches the python script trackGUI.py
os.system('start "Controls trackGui.py" python trackGUI.py 3')
