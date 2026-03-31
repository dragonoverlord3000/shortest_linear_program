import json

with open("../build/benchmarks/full.json", "r") as f:
    j_bp = json.load(f)

with open("../build/benchmarks/full_gp.json", "r") as f:
    j_gp = json.load(f)


scheme2add_bp = [
    {
        "scheme_name": scheme["scheme_name"],
        "U": scheme["results"]["U"]["best_add"],
        "V": scheme["results"]["V"]["best_add"],
        "W": scheme["results"]["W"]["best_add"],
    }
    for scheme in j_bp["results"][0]["details"]["schemes"]
]

scheme2add_gp = [
    {
        "scheme_name": scheme["scheme_name"],
        "U": scheme["results"]["U"]["best_add"],
        "V": scheme["results"]["V"]["best_add"],
        "W": scheme["results"]["W"]["best_add"],
    }
    for scheme in j_gp["results"][0]["details"]["schemes"]
]

best_scheme = {}
best_val = float("inf")

print("combining like schemes")
schemes = []
for scheme_bp in scheme2add_bp:
    for scheme_gp in scheme2add_gp:
        if scheme_bp["scheme_name"] != scheme_gp["scheme_name"]:
            continue
        schemes.append((scheme_bp, scheme_gp))

cost2cnt = {}
print(f"going through {len(scheme2add_bp)} schemes")
for scheme_bp, scheme_gp in schemes:
    cur = (
        min(scheme_bp["U"], scheme_gp["U"])
        + min(scheme_bp["V"], scheme_gp["V"])
        + min(scheme_bp["W"], scheme_gp["W"])
    )
    cost2cnt[cur] = cost2cnt.get(cur, 0) + 1

    if cur < best_val:
        best_scheme = (scheme_bp, scheme_gp)
        best_val = cur

print("cost2cnt: ", cost2cnt)
print("best scheme: ", best_scheme)
