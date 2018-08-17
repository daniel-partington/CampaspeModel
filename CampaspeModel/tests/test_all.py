import os
import pkgutil
import subprocess
import sys


def run_model(folder, name):
    os.chdir(folder)
    return subprocess.check_output([pyexe, name])


pyexe = sys.executable
if "pythonw" in pyexe:
    pyexe = pyexe.replace("pythonw", "python")

this_dir = os.getcwd()
model_dir = os.path.join(this_dir, "..", "models")
os.chdir(model_dir)


# Get all CampaspeModels
model_pkgs = [name for _, name, _ in pkgutil.walk_packages([model_dir])]
# Shorten model package list because I don't have the other models and cannot build them
model_pkgs = ['GW_link_Integrated']


def test_GW_link_Integrated_run_forecast():
    expected_lines = [
        "swgw_exchanges",
        "avg_depth_to_gw",
        "ecol_depth_to_gw",
        "trigger_heads"
    ]

    folder = 'GW_link_Integrated'  # TODO: Replace this with path fixture

    # ideally we'd run the stub and return the model object to pull the results from the output directly
    # from CampaspeModel.models.GW_link_Integrated import GW_link_Integrated_run_forecast as forecast
    # modflow_model = forecast.main()
    results = run_model(folder, 'GW_link_Integrated_run_forecast.py')

    for line in expected_lines:
        assert line in results
    # End for

# End test_GW_link_Integrated_run_forecast()
