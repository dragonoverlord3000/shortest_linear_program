import sympy as sp
from collections import defaultdict, deque

basis_change = sp.Matrix(
    [
        [0, 0, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 1, 0],
        [0, 0, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0, -1],
        [-1, 0, 0, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, -1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 1],
        [1, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 1, 0, 0, 0, 0],
    ]
)

scheme = """
s_40 w_9 -1
s_23 s_40 1
s_37 s_40 -1
s_16 s_37 1
s_33 s_37 -1
s_27 s_33 -1
s_6 s_27 1
s_19 s_27 -1
s_29 s_33 -1
s_9 s_29 1
s_24 s_29 -1
s_3 s_24 -1
s_10 s_24 -1
s_46 w_7 -1
s_33 s_46 1
s_7 s_46 -1
s_52 w_6 -1
s_39 s_52 -1
s_22 s_39 1
s_36 s_39 -1
s_26 s_36 -1
s_14 s_26 -1
s_18 s_26 -1
s_32 s_36 -1
s_23 s_32 1
s_27 s_32 -1
s_21 s_52 -1
s_43 w_5 -1
s_20 s_43 -1
s_41 s_43 -1
s_5 s_41 -1
s_38 s_41 -1
s_2 s_38 1
s_35 s_38 -1
s_1 s_35 1
s_31 s_35 -1
s_13 s_31 1
s_26 s_31 -1
s_45 w_4 -1
s_44 s_45 -1
s_19 s_44 1
s_42 s_44 -1
s_31 s_42 1
s_41 s_42 1
s_4 s_45 -1
s_49 w_3 -1
s_34 s_49 -1
s_18 s_34 1
s_30 s_34 -1
s_16 s_30 -1
s_24 s_30 1
s_12 s_49 -1
s_48 w_2 -1
s_38 s_48 -1
s_11 s_48 -1
s_50 w_8 -1
s_47 s_50 -1
s_28 s_47 -1
s_20 s_28 1
s_22 s_28 -1
s_8 s_47 1
s_15 s_50 -1
s_51 w_1 -1
s_25 s_51 -1
s_1 s_25 1
s_10 s_25 1
s_17 s_51 -1
"""

# GET ORDERING START

depends_on = {}
inv_depends_on = {}
for line in scheme.split("\n"):
    if len(line) < 3:
        continue
    v1, v2, _ = line.split(" ")
    if v1 not in depends_on:
        depends_on[v1] = set()
        inv_depends_on[v1] = set()
    if v2 not in depends_on:
        depends_on[v2] = set()
        inv_depends_on[v2] = set()

    depends_on[v2].add(v1)
    inv_depends_on[v1].add(v2)

ordering = []
in_degree = {node: len(values) for node, values in depends_on.items()}
queue = deque()
for node in in_degree:
    if not in_degree[node]:
        queue.append(node)

print(in_degree)

while len(queue):
    node = queue.popleft()
    if node not in ordering:
        ordering.append(node)

    for neigh in inv_depends_on[node]:
        in_degree[neigh] -= 1
        if not in_degree[neigh]:
            queue.append(neigh)

# GET ORDERING END

T = "u"
if "v" in scheme:
    T = "v"
elif "w" in scheme:
    T = "w"

print(f"T = {T}")

num_inputs = 23 if T == "w" else 9
name = {"u": "\\tilde{A}", "v": "\\tilde{B}", "w": "M"}[T]

var2sum = {}
for line in scheme.split("\n"):
    if len(line) < 3:
        continue
    v1, v2, sgn = line.split(" ")
    var2sum[v2] = var2sum.get(v2, []) + [(sgn, v1)]

var_map = {}
for idx in range(num_inputs):
    if T == "w":
        var_map[f"s_{idx + 1}"] = "M_{" + f"{idx + 1}" + "}"
    else:
        i = idx // 3 + 1
        j = idx % 3 + 1
        var_map[f"s_{idx + 1}"] = f"{name}_" + "{" + f"{i},{j}" + "}"

seen = set(var_map.keys())
not_seen = set()
for line in scheme.split("\n"):
    if len(line) < 3:
        continue
    v1, v2, _ = line.split(" ")
    for v in [v1, v2]:
        if "s" == v[0] and v not in seen:
            not_seen.add(v)
        if v[0] in "uvw" and v not in var_map:
            var_map[v] = v.split("_")[0] + "_{" + v.split("_")[1] + "}"

not_seen_ordered = []
for node in ordering:
    if node in not_seen:
        not_seen_ordered.append(node)

print(ordering)
print(not_seen)
print(not_seen_ordered)

cnt = 0
for v in not_seen_ordered:
    cnt += 1
    var_map[v] = "t_{" + f"{cnt}" + "}"
    seen.add(v)

for k, v in var_map.items():
    if "t_" in v:
        var_map[k] = v.replace("t", "s")

equations = []
for og_var, S in var2sum.items():
    var = var_map.get(og_var, og_var)
    ans = var
    ans += " = "

    for sgn, v in S:
        v = var_map.get(v, v)
        if int(sgn) == 1:
            if ans[-3:] == " = ":
                ans += f" {v}"
            else:
                ans += f" + {v}"
        else:
            if ans[-3:] == " = ":
                ans += f"-{v}"
            else:
                ans += f" - {v}"
    equations.append((ans[0], int(ans.split("{")[1].split("}")[0]), ans))
equations.sort(key=lambda x: (x[0], x[1]))
print(equations)

for i in range(len(equations)):
    _, _, ans = equations[i]
    print(f"{ans} \\\\")
