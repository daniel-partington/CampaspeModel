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
        "('swgw_exchanges', rec.array([( 866104.125,  688921.5,  1748.95361328)],",
        "dtype=[('406201', '<f8'), ('406202', '<f8'), ('406265', '<f8')]))",
        "('avg_depth_to_gw', rec.array([ ( 18.57427216,  8.39050674,  7.42168935,  10.43254089,  8.08312225,  3.61949921,  18.96347264,  11.81394644,  11.36267471,  11.83951416,  10.85602529,  7.32978058)],",
        "dtype=[('1', '<f8'), ('2', '<f8'), ('3', '<f8'), ('4', '<f8'), ('5', '<f8'), ('6', '<f8'), ('7', '<f8'), ('8', '<f8'), ('9', '<f8'), ('10', '<f8'), ('11', '<f8'), ('12', '<f8')]))",
        "('ecol_depth_to_gw', rec.array([( 11.37800598,  7.33989716,  18.57427216)],",
        "dtype=[('62599', '<f8'), ('47249', '<f8'), ('79329', '<f8')]))",
        "('trigger_heads', rec.array([( 121.24111176,  85.66869354)],",
        "dtype=[('62589', '<f8'), ('79324', '<f8')]))",
        "('swgw_exchanges', rec.array([( 848939.5,  690552.4375,  1748.95361328)],",
        "dtype=[('406201', '<f8'), ('406202', '<f8'), ('406265', '<f8')]))",
        "('avg_depth_to_gw', rec.array([ ( 18.5730896,  8.39127477,  7.46096802,  10.38856888,  8.08272552,  3.63551712,  18.96313259,  11.81394869,  11.36324183,  11.82771835,  10.83888204,  7.33062363)],",
        "dtype=[('1', '<f8'), ('2', '<f8'), ('3', '<f8'), ('4', '<f8'), ('5', '<f8'), ('6', '<f8'), ('7', '<f8'), ('8', '<f8'), ('9', '<f8'), ('10', '<f8'), ('11', '<f8'), ('12', '<f8')]))",
        "('ecol_depth_to_gw', rec.array([( 11.34143066,  7.34079742,  18.5730896)],",
        "dtype=[('62599', '<f8'), ('47249', '<f8'), ('79329', '<f8')]))",
        "('trigger_heads', rec.array([( 121.24082184,  85.66912079)],",
        "dtype=[('62589', '<f8'), ('79324', '<f8')]))"
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
