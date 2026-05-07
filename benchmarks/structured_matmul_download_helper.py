import urllib.request
import sympy as sp
import json
import os

from sympy.parsing.sympy_parser import (
    parse_expr,
    standard_transformations,
    implicit_multiplication_application,
)

transformations = standard_transformations + (implicit_multiplication_application,)


# List generated with AI
# the files from https://github.com/khoruzhii/flip-cpd/tree/main/data/schemes_paper
txt_files = [
    "gg-222-rank7-rec-0-0-0-z.txt",
    "gg-333-rank23-rec-0-0-0-z.txt",
    "gg-444-rank49-rec-0-0-0-z.txt",
    "gg-555-rank97-rec-0-0-0-z.txt",
    "gt-333-rank17-rec-0-0-0-z.txt",
    "gt-333-rank17-rec-3-0-0-q.txt",
    "gt-333-rank17-rec-9-0-0-z.txt",
    "gt-444-rank34-rec-1-0-0-q.txt",
    "gt-444-rank34-rec-4-0-0-z.txt",
    "gt-444c-rank26-rec-2-0-0-z-4.txt",
    "gt-555-rank63-rec-1-0-0-z.txt",
    "gt-555-rank74-rec-6-0-0-q.txt",
    "kg-333-rank14-rec-0-0-0-q.txt",
    "kg-333-rank15-rec-0-0-0-z.txt",
    "kg-444-rank36-rec-0-0-0-q.txt",
    "kg-444-rank36-rec-0-0-0-z.txt",
    "kg-555-rank70-rec-0-0-0-z.txt",
    "kk-222-rank1-rec-0-0-0-z.txt",
    "kk-333-rank9-rec-0-0-0-z.txt",
    "kk-444-rank20-rec-0-0-0-q.txt",
    "kk-444-rank22-rec-0-0-0-z.txt",
    "kk-555-rank48-rec-0-0-0-z.txt",
    "kt-222-rank1-rec-0-0-0-z.txt",
    "kt-333-rank6-rec-0-0-0-z.txt",
    "kt-444-rank15-rec-0-0-0-z.txt",
    "kt-555-rank30-rec-0-0-0-z.txt",
    "sg-222-rank6-rec-0-0-0-z.txt",
    "sg-222-rank6-rec-1-0-0-q.txt",
    "sg-333-rank18-rec-0-0-0-z.txt",
    "sg-444-rank40-rec-0-0-0-z.txt",
    "sg-555-rank75-rec-0-0-0-z.txt",
    "sk-222-rank3-rec-0-2-0-z.txt",
    "sk-333-rank11-rec-0-0-0-z.txt",
    "sk-333-rank11-rec-0-2-0-q.txt",
    "sk-444-rank24-rec-0-0-0-q.txt",
    "sk-444-rank26-rec-0-0-0-z.txt",
    "sk-555-rank50-rec-0-0-0-z.txt",
    "ss-222-rank5-rec-0-2-0-q.txt",
    "ss-222-rank6-rec-0-5-0-z.txt",
    "ss-222-rank6-rec-2-2-0-z.txt",
    "ss-333-rank14-rec-0-2-0-q.txt",
    "ss-333-rank15-rec-0-0-0-z.txt",
    "ss-444-rank30-rec-0-1-0-q.txt",
    "ss-444-rank32-rec-0-0-0-z.txt",
    "ss-444-rank33-rec-0-0-0-z.txt",
    "ss-555-rank59-rec-0-0-0-z.txt",
    "st-222-rank4-rec-0-2-0-z.txt",
    "st-333-rank10-rec-0-3-0-q.txt",
    "st-333-rank11-rec-0-2-1-z.txt",
    "st-444-rank20-rec-0-4-0-q.txt",
    "st-444-rank22-rec-0-2-2-z.txt",
    "st-555-rank35-rec-0-3-0-q.txt",
    "st-555-rank39-rec-0-2-1-z.txt",
    "sw-222-rank5-rec-0-1-0-q.txt",
    "sw-222-rank6-rec-0-2-3-z.txt",
    "sw-222-rank6-rec-0-3-2-z.txt",
    "sw-222-rank6-rec-2-1-1-z.txt",
    "sw-333-rank14-rec-0-2-1-q.txt",
    "sw-333-rank15-rec-0-0-0-z.txt",
    "sw-444-rank31-rec-0-0-0-q.txt",
    "sw-444-rank35-rec-0-0-1-z.txt",
    "sw-444-rank35-rec-0-1-0-z.txt",
    "sw-555-rank64-rec-0-0-0-z.txt",
    "sw-555-rank64-rec-0-2-0-q.txt",
    "ug-333-rank17-rec-6-0-0-z.txt",
    "ug-333-rank17-rec-8-0-0-q.txt",
    "ug-333-rank17-rec-8-0-0-z.txt",
    "ug-444-rank34-rec-12-0-0-q.txt",
    "ug-444-rank34-rec-12-0-0-z.txt",
    "ug-555-rank63-rec-17-0-0-z.txt",
    "uk-333-rank10-rec-0-3-0-q.txt",
    "uk-333-rank11-rec-0-5-0-z.txt",
    "uk-444-rank24-rec-0-6-0-q.txt",
    "uk-444-rank24-rec-0-6-0-z.txt",
    "uk-555-rank45-rec-0-7-0-z.txt",
    "uk-555-rank46-rec-0-9-0-z.txt",
    "ul-333-rank13-rec-0-2-5-z.txt",
    "ul-333-rank13-rec-0-3-5-q.txt",
    "ul-333-rank13-rec-0-3-5-z.txt",
    "ul-333-rank13-rec-0-4-4-z.txt",
    "ul-333-rank13-rec-0-5-3-z.txt",
    "ul-333-rank13-rec-1-3-4-z.txt",
    "ul-333-rank13-rec-1-4-3-z.txt",
    "ul-333-rank13-rec-2-2-3-z.txt",
    "ul-333-rank13-rec-2-3-2-z.txt",
    "ul-333-rank13-rec-2-3-3-q.txt",
    "ul-333-rank13-rec-2-3-3-z.txt",
    "ul-333-rank13-rec-3-2-2-q.txt",
    "ul-333-rank13-rec-3-2-2-z.txt",
    "ul-444-rank27-rec-0-5-8-q.txt",
    "ul-444-rank27-rec-0-5-8-z.txt",
    "ul-444-rank27-rec-0-6-7-q.txt",
    "ul-444-rank27-rec-0-6-7-z.txt",
    "ul-444-rank27-rec-0-7-6-q.txt",
    "ul-444-rank27-rec-0-7-6-z.txt",
    "ul-444-rank27-rec-0-8-5-q.txt",
    "ul-444-rank27-rec-0-8-5-z.txt",
    "ul-444-rank27-rec-1-4-7-q.txt",
    "ul-444-rank27-rec-1-5-7-z.txt",
    "ul-444-rank27-rec-1-6-6-q.txt",
    "ul-444-rank27-rec-1-6-6-z.txt",
    "ul-444-rank27-rec-1-7-5-q.txt",
    "ul-444-rank27-rec-1-7-5-z.txt",
    "ul-444-rank27-rec-2-4-6-q.txt",
    "ul-444-rank27-rec-2-5-5-q.txt",
    "ul-444-rank27-rec-2-5-6-z.txt",
    "ul-444-rank27-rec-2-6-2-q.txt",
    "ul-444-rank27-rec-2-6-5-z.txt",
    "ul-444-rank27-rec-3-4-5-q.txt",
    "ul-444-rank27-rec-3-4-5-z.txt",
    "ul-444-rank27-rec-3-5-4-z.txt",
    "ul-444-rank27-rec-4-3-4-z.txt",
    "ul-555-rank47-rec-0-10-9-z.txt",
    "ul-555-rank47-rec-0-11-8-z.txt",
    "ul-555-rank47-rec-0-8-11-z.txt",
    "ul-555-rank47-rec-0-9-10-z.txt",
    "ul-555-rank47-rec-1-10-7-z.txt",
    "ul-555-rank47-rec-1-8-10-z.txt",
    "ul-555-rank47-rec-1-9-9-q.txt",
    "ul-555-rank47-rec-2-6-9-z.txt",
    "ul-555-rank47-rec-2-7-9-q.txt",
    "ul-555-rank47-rec-2-9-8-z.txt",
    "ul-555-rank47-rec-3-5-7-z.txt",
    "ul-555-rank47-rec-3-7-6-z.txt",
    "ul-555-rank47-rec-3-7-7-q.txt",
    "us-222-rank5-rec-0-3-0-z.txt",
    "us-222-rank5-rec-1-2-0-z.txt",
    "us-222-rank5-rec-2-1-0-q.txt",
    "us-333-rank14-rec-0-6-0-z.txt",
    "us-333-rank14-rec-1-3-0-z.txt",
    "us-333-rank14-rec-1-4-1-q.txt",
    "us-333-rank14-rec-1-5-0-q.txt",
    "us-333-rank14-rec-1-5-0-z.txt",
    "us-333-rank14-rec-2-3-0-q.txt",
    "us-444-rank28-rec-0-7-0-q.txt",
    "us-444-rank28-rec-1-6-0-q.txt",
    "us-444-rank29-rec-0-9-0-z.txt",
    "us-444-rank29-rec-1-8-0-z.txt",
    "us-555-rank52-rec-0-12-0-z.txt",
    "us-555-rank52-rec-1-11-0-z.txt",
    "ut-444-rank19-rec-0-8-1-z.txt",
    "ut-444-rank19-rec-0-8-2-q.txt",
    "ut-444-rank19-rec-0-9-1-q.txt",
    "ut-444-rank19-rec-0-9-3-z.txt",
    "ut-444-rank19-rec-1-7-2-q.txt",
    "ut-444-rank19-rec-1-7-3-z.txt",
    "ut-444-rank19-rec-1-8-1-q.txt",
    "ut-444-rank19-rec-1-8-2-z.txt",
    "ut-444-rank19-rec-2-6-1-q.txt",
    "ut-444-rank19-rec-2-7-2-z.txt",
    "ut-444-rank19-rec-3-6-3-z.txt",
    "ut-444-rank19-rec-4-4-2-z.txt",
    "ut-444c-rank14-rec-2-4-3-z-5.txt",
    "ut-555-rank32-rec-0-10-1-q.txt",
    "ut-555-rank32-rec-0-10-3-z.txt",
    "ut-555-rank32-rec-0-11-0-q.txt",
    "ut-555-rank32-rec-0-11-2-z.txt",
    "ut-555-rank32-rec-1-10-0-q.txt",
    "ut-555-rank32-rec-1-10-2-z.txt",
    "ut-555-rank32-rec-1-7-4-z.txt",
    "ut-555-rank32-rec-2-7-2-q.txt",
    "ut-555-rank32-rec-2-8-1-q.txt",
    "ut-555-rank32-rec-2-8-3-z.txt",
    "ut-555-rank32-rec-2-9-1-z.txt",
    "ut-555-rank32-rec-3-8-2-z.txt",
    "uu-444-rank19-rec-0-14-0-q.txt",
    "uu-444-rank19-rec-0-14-0-z.txt",
    "uu-444-rank19-rec-1-12-0-z.txt",
    "uu-444-rank19-rec-1-13-0-q.txt",
    "uu-444-rank19-rec-1-13-0-z.txt",
    "uu-444-rank19-rec-2-12-0-q.txt",
    "uu-444-rank19-rec-2-12-0-z.txt",
    "uu-444-rank19-rec-2-9-0-z.txt",
    "uu-444-rank19-rec-3-11-0-q.txt",
    "uu-444-rank19-rec-3-11-0-z.txt",
    "uu-444-rank19-rec-4-10-0-q.txt",
    "uu-444-rank19-rec-4-10-0-z.txt",
    "uu-555-rank32-rec-1-16-0-q.txt",
    "uu-555-rank32-rec-1-16-0-z.txt",
    "uu-555-rank32-rec-2-15-0-q.txt",
    "uu-555-rank32-rec-2-15-0-z.txt",
    "uu-555-rank32-rec-3-14-0-q.txt",
    "uu-555-rank32-rec-3-14-0-z.txt",
    "uu-555-rank32-rec-4-13-0-q.txt",
    "uu-555-rank32-rec-4-13-0-z.txt",
    "uu-555-rank32-rec-5-10-0-q.txt",
    "uu-555-rank32-rec-5-12-0-z.txt",
    "uw-222-rank5-rec-0-3-0-z.txt",
    "uw-222-rank5-rec-1-2-0-z.txt",
    "uw-222-rank5-rec-2-1-0-q.txt",
    "uw-333-rank14-rec-0-2-1-q.txt",
    "uw-333-rank14-rec-0-6-0-z.txt",
    "uw-333-rank14-rec-1-5-0-z.txt",
    "uw-333-rank14-rec-2-3-0-q.txt",
    "uw-444-rank29-rec-0-9-0-q.txt",
    "uw-444-rank29-rec-0-9-0-z.txt",
    "uw-444-rank29-rec-1-8-0-q.txt",
    "uw-444-rank29-rec-1-8-0-z.txt",
    "uw-444-rank29-rec-2-6-0-q.txt",
    "uw-555-rank54-rec-0-13-0-q.txt",
    "uw-555-rank54-rec-0-9-0-z.txt",
    "uw-555-rank54-rec-1-12-0-q.txt",
    "wg-222-rank6-rec-0-0-0-z.txt",
    "wg-222-rank6-rec-1-0-0-q.txt",
    "wg-333-rank18-rec-0-0-0-z.txt",
    "wg-333-rank18-rec-2-0-0-q.txt",
    "wg-444-rank40-rec-0-0-0-q.txt",
    "wg-444-rank40-rec-0-0-0-z.txt",
    "wg-555-rank75-rec-0-0-0-z.txt",
    "wt-222-rank4-rec-0-2-0-z.txt",
    "wt-333-rank11-rec-0-2-1-z.txt",
    "wt-333-rank11-rec-0-3-1-q.txt",
    "wt-444-rank22-rec-0-2-0-z.txt",
    "wt-444-rank22-rec-0-2-2-q.txt",
    "wt-555-rank40-rec-0-1-2-q.txt",
    "wt-555-rank40-rec-0-2-1-q.txt",
    "wt-555-rank44-rec-0-4-2-z.txt",
    "ww-222-rank5-rec-0-2-0-q.txt",
    "ww-222-rank6-rec-0-5-0-z.txt",
    "ww-222-rank6-rec-2-2-0-z.txt",
    "ww-333-rank15-rec-0-0-0-z.txt",
    "ww-333-rank15-rec-0-3-0-q.txt",
    "ww-333-rank15-rec-3-0-0-q.txt",
    "ww-444-rank33-rec-0-0-0-q.txt",
    "ww-444-rank35-rec-0-1-0-z.txt",
    "ww-555-rank68-rec-0-1-0-z.txt",
    "ww-555-rank68-rec-0-2-0-q.txt",
]

os.makedirs("benchmarks/structured_matmul", exist_ok=True)
output_dir = "benchmarks/structured_matmul/dataset"
os.makedirs(output_dir, exist_ok=True)

url_base = "https://raw.githubusercontent.com/khoruzhii/flip-cpd/refs/heads/main/data/schemes_paper/"
for file_name in txt_files:
    download_url = url_base + file_name
    destination = os.path.join(output_dir, file_name)
    print(f"Downloading {file_name}...")
    urllib.request.urlretrieve(download_url, destination)
    print(f"{file_name} downloaded")

    scheme_name = file_name[:-4]
    r = int(file_name.split("rank")[1].split("-")[0])
    n = int(file_name[3])
    type_a = file_name[0]
    type_b = file_name[1]

    A = [[0] * n for _ in range(n)]
    B = [[0] * n for _ in range(n)]

    def set_matrix(M, type_m, name):
        at = 1
        if type_m == "g":
            for i in range(n):
                for j in range(n):
                    M[i][j] = sp.Symbol(f"{name}{at}")
                    at += 1
        elif type_m == "s":
            for i in range(n):
                for j in range(i, n):
                    M[j][i] = M[i][j] = sp.Symbol(f"{name}{at}")
                    at += 1
        elif type_m == "k":
            for i in range(n):
                for j in range(i + 1, n):
                    M[i][j] = sp.Symbol(f"{name}{at}")
                    M[j][i] = -M[i][j]
                    at += 1
        elif type_m == "u":
            for i in range(n):
                for j in range(i, n):
                    M[i][j] = sp.Symbol(f"{name}{at}")
                    at += 1
        elif type_m == "l":
            for j in range(n):
                for i in range(j, n):
                    M[i][j] = sp.Symbol(f"{name}{at}")
                    at += 1
        elif type_m == "w":
            for i in range(n):
                for j in range(i, n):
                    x = sp.Symbol(f"{name}{at}")
                    M[i][j] = x
                    if i != j:
                        M[j][i] = -x
                    at += 1
        return M


    A = set_matrix(A, type_a, "a")
    if type_b != "t":
        B = set_matrix(B, type_b, "b")
        B = sp.Matrix(B)
    else:
        B = set_matrix(B, type_a, "b")
        B = sp.Matrix(B).T

    A = sp.Matrix(A)


    def num_vars(type_m, n):
        if type_m == "g":
            return n * n
        if type_m in {"s", "u", "l", "w"}:
            return n * (n + 1) // 2
        if type_m == "k":
            return n * (n - 1) // 2
        raise ValueError(f"Unknown type: {type_m}")

    num_a = num_vars(type_a, n)

    if type_b == "t":
        num_b = num_a
    else:
        num_b = num_vars(type_b, n)

    B2A = {sp.Symbol(f"b{at}"): sp.Symbol(f"a{at}") for at in range(1, num_b + 1)}

    var2num_a = {sp.Symbol(f"a{v}"): v for v in range(1, num_a + 1)}
    var2num_b = {sp.Symbol(f"b{v}"): v for v in range(1, num_b + 1)}
    var2num_m = {sp.Symbol(f"m{v}"): v for v in range(1, r + 1)}


    local_dict = {f"m{v}": sp.Symbol(f"m{v}") for v in range(1, r + 1)}
    for i in range(1, n**2 + 1):
        local_dict[f"a{i}"] = sp.Symbol(f"a{i}")
        local_dict[f"b{i}"] = sp.Symbol(f"b{i}")

    with open(f"{output_dir}/{file_name}", "r") as f:
        data = f.read()

    eqs_a = []
    eqs_b = []
    eqs_m = []
    m_flag = False
    for line in data.split("\n"):
        if not line:
            m_flag = True
            continue

        if m_flag:
            lhs, rhs = line.split("=")
            rhs = rhs.strip(" ")
            rhs = parse_expr(rhs, local_dict=local_dict, transformations=transformations)
            eqs_m.append((line, rhs))

        else:
            lhs, rhs = line.split("=")
            eq1, eq2 = rhs.split(")(")
            eq1 = eq1.replace("(", "").replace(")", "").strip(" ")
            eq2 = eq2.replace("(", "").replace(")", "").strip(" ")
            eq1 = parse_expr(eq1, local_dict=local_dict, transformations=transformations)
            eq2 = parse_expr(eq2, local_dict=local_dict, transformations=transformations)

            eqs_a.append((line, eq1))
            if type_b == "t":
                eqs_a.append((line, eq2.subs(B2A)))
            else:
                eqs_b.append((line, eq2))


    def make_binary_matrix(eqs, variable2num):
        mat = []
        for _, eq in eqs:
            row = [0] * len(variable2num)
            for var, num in variable2num.items():
                row[num-1] = eq.coeff(var) % 2
            mat.append(row)
        return mat

    # info dump
    with open(f"{output_dir}/{scheme_name}.json", "w") as f:
        json.dump(
            {
                "scheme name": scheme_name,
                "A": [[str(A[i, j]) for j in range(n)] for i in range(n)],
                "B": [[str(B[i, j]) for j in range(n)] for i in range(n)],
                "type_a": type_a,
                "type_b": type_b,
                "original_scheme": data,
            }, f, indent=2
        )

    # the actual instances
    with open(f"{output_dir}/{scheme_name}_A.txt", "w") as f:
        mat = make_binary_matrix(eqs_a, var2num_a)
        mat_a_s = f"{len(mat)} {len(mat[0])}"
        for row in mat:
            mat_a_s += "\n"
            for j in range(len(row)):
                mat_a_s += str(row[j])
                if j < len(row) - 1:
                    mat_a_s += " "
        f.write(mat_a_s)

    if eqs_b:
        with open(f"{output_dir}/{scheme_name}_B.txt", "w") as f:
            mat = make_binary_matrix(eqs_b, var2num_b)
            mat_b_s = f"{len(mat)} {len(mat[0])}"
            for row in mat:
                mat_b_s += "\n"
                for j in range(len(row)):
                    mat_b_s += str(row[j])
                    if j < len(row) - 1:
                        mat_b_s += " "
            f.write(mat_b_s)


    with open(f"{output_dir}/{scheme_name}_M.txt", "w") as f:
        mat = make_binary_matrix(eqs_m, var2num_m)
        mat_m_s = f"{len(mat)} {len(mat[0])}"
        for row in mat:
            mat_m_s += "\n"
            for j in range(len(row)):
                mat_m_s += str(row[j])
                if j < len(row) - 1:
                    mat_m_s += " "
        f.write(mat_m_s)

print("All files downloaded and converted successfully!")

