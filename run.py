import os

graphs=["NY_10","NY_15","NY_20","NY_25",
        "BAY_10","BAY_15","BAY_20","BAY_25",
        "COL_10","COL_15","COL_20","COL_25",
        "FLA_10","FLA_15","FLA_20","FLA_25"]

types=["V0","V8","V11","naive"]
thresholds = ["0"]
result_path = "./result_road_variable_w/"

for graph in graphs:
    for type_ in types:
        for threshold in thresholds:
            command = "/usr/bin/time -v timeout 20h ./wcsd_626 road_networks/"+ graph + "_tree txt road " + type_ + " " + threshold + " > " + result_path + graph + "_tree_" + type_ + "_" + threshold + " 2>&1"
            print(command)
            os.system(command)