import os

graphs=["eu-2005_w15","eu-2005_w20",
        "movielens_25m_w15","movielens_25m_w20",
        "uk-2007_w15","uk-2007_w20"]

types=["V0","V8","naive"]
thresholds = ["10"]
result_path = "./result_social_variable_w/"

for graph in graphs:
    for type_ in types:
        for threshold in thresholds:
            command = "/usr/bin/time -v timeout 30h ./wcsd_626 "+ graph + "_degree bin social " + type_ + " " + threshold + " > " + result_path + graph + "_hybrid_" + type_ + "_" + threshold + " 2>&1"
            print(command)
            os.system(command)