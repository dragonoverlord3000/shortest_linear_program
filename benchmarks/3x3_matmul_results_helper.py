import json
import os

adds2print = {60}

scheme2adds = []
scheme_names = set()

for json_file in os.listdir("../build/benchmarks/"):
    if json_file[-5:] != ".json" or not "3x3" in json_file:
        continue
    with open(f"../build/benchmarks/{json_file}", "r") as f:
        j = json.load(f)
        scheme2add = {
            scheme["scheme_name"]: {
                "U": scheme["results"]["U"]["best_add"]
                if "U" in scheme["results"]
                else float("inf"),
                "V": scheme["results"]["V"]["best_add"]
                if "V" in scheme["results"]
                else float("inf"),
                "W": scheme["results"]["W"]["best_add"]
                if "W" in scheme["results"]
                else float("inf"),
                "search_method": j["config"]["search_method"],
            }
            for scheme in j["results"][0]["details"]["schemes"]
        }
        for scheme in j["results"][0]["details"]["schemes"]:
            scheme_names.add(scheme["scheme_name"])
        scheme2adds.append(scheme2add)

schemes = {}
best_val = float("inf")

print("combining like schemes")
for scheme_name in scheme_names:
    for scheme2add in scheme2adds:
        scheme = scheme2add[scheme_name]
        if scheme_name not in schemes:
            schemes[scheme_name] = scheme
            schemes[scheme_name]["search_method"] = {
                schemes[scheme_name]["search_method"]
            }
            continue

        for type in ["U", "V", "W"]:
            if scheme[type] < schemes[scheme_name][type]:
                schemes[scheme_name][type] = scheme[type]
                schemes[scheme_name]["search_method"].add(scheme["search_method"])

cost2cnt = {}
for scheme in schemes.values():
    cost = 0
    for type in ["U", "V", "W"]:
        cost += scheme[type]
    cost2cnt[cost] = cost2cnt.get(cost, 0) + 1

print("cost2cnt: ", cost2cnt)
print("best cost: ", min(cost2cnt.keys()))

for scheme_name, scheme in schemes.items():
    cost = 0
    for type in ["U", "V", "W"]:
        cost += scheme[type]
    if cost in adds2print:
        print(f"{scheme_name}: {cost}")
