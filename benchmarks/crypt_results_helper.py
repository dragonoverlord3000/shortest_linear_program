import json
import os

scheme2add = {}
scheme_names = set()

for json_file in os.listdir("../build/benchmarks/"):
    if json_file[-5:] != ".json" or not "crypt" in json_file:
        continue
    with open(f"../build/benchmarks/{json_file}", "r") as f:
        j = json.load(f)
        for scheme in j["results"][0]["details"]["schemes"]:
            scheme2add[scheme["scheme_name"]] = min(
                scheme2add.get(scheme["scheme_name"], float("inf")),
                scheme["results"]["best_add"],
            )

print(scheme2add)
