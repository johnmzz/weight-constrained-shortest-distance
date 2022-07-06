import random

def parse_road(path, name, delim, num_weights):
    sizes = {"movielens_10m": (80555,10000054),
             "eu-2005": (862664, 16138468),
             "eswiki-2013": (970331, 21184931),
             "movielens_25m": (221588, 25000095),
             "frwiki": (1350986, 31037302),
             "uk-2007": (1000000, 37061970),
    }

    f = open("./data/social_networks_org/" + name + ".txt", 'r')
    fgraph = open("./data/social_networks_org/" + name + "_" + str(num_weights) + ".txt", 'w+')

    fgraph.write(str(sizes[name][0]) + " " + str(sizes[name][1]) + " " + str(num_weights) + "\n")
    for line in f.readlines():
        lst = line.split(delim)

        u = int(lst[0])
        v = int(lst[1])

        fgraph.write(str(u) + " " + str(v) + " " + str(random.randint(1,num_weights)) + "\n")


graphs = ["movielens_10m","eu-2005","eswiki-2031","movielens_25m","frwiki","uk-2007"]
path = "./data/social_networks_org/"

for graph in graphs:
    print(f"parsing {graph}...")
    parse_road(path, graph, " ", 10)