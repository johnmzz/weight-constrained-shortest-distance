import random

def parse_road(path, name, delim, num_weights):
    sizes = {"NY": (264346, 733846),
             "BAY": (321270, 800172),
             "COL": (435666, 1057066),
             "FLA": (1070376, 2712798),
             "CAL": (1890815, 4657742),
             "E": (3598623, 8778114),
             "W": (6262104, 15248146),
             "CTR": (14081816, 34292496),
             "USA": (23947347, 58333344)
    }

    f = open("./data/road_networks_org/" + name + ".txt", 'r')
    fgraph = open("./data/road_networks/" + name + "_" + str(num_weights) + ".txt", 'w+')

    fgraph.write(str(sizes[name][0]) + " " + str(sizes[name][1]) + " " + str(num_weights) + "\n")
    for line in f.readlines():
        lst = line.split(delim)

        u = int(lst[1]) - 1
        v = int(lst[2]) - 1

        fgraph.write(str(u) + " " + str(v) + " " + str(random.randint(1,num_weights)) + "\n")


graphs = ["NY","BAY","COL","FLA","CAL","E","W","CTR","USA"]
path = "./data/road_networks"

for graph in graphs:
    print(f"parsing {graph}...")
    parse_road(path, graph, " ", 25)