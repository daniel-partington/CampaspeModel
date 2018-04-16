import sys
import subprocess
import os

orig = os.getcwd()

pyexe = sys.executable
if "pythonw" in pyexe:
    pyexe = pyexe.replace("pythonw", "python")

from steady_state import Campaspe_SS_build
from transient import Campaspe_transient_build

os.chdir(orig)
def run_model(folder, name):
    os.chdir(folder)
    print(subprocess.check_output([pyexe, name]))

models = [("./steady_state", "Campaspe_SS_run_flow.py"), 
          ("./steady_state", "Campaspe_SS_run_transport.py"),
          ("./transient", "Campaspe_transient_run_flow.py"),
          ("./transient", "Campaspe_transient_run_transport_C14.py"),
          ("./transient", "Campaspe_transient_run_transport_RnEC.py"),]     

for model in models:
    os.chdir(orig)
    run_model(*model)

