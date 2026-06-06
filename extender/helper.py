import sympy as sp

flip = True
B = [
    [-1, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, -1, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, -1, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, -1, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, -1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, -1, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, -1, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, -1, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, -1],
]

# U = 9, V = 9, W = 23
num_inputs = 23

scheme = """
s_40 w_9 -1
s_23 s_40 -1
s_37 s_40 -1
s_16 s_37 -1
s_33 s_37 -1
s_27 s_33 -1
s_6 s_27 -1
s_19 s_27 1
s_29 s_33 -1
s_9 s_29 -1
s_24 s_29 -1
s_3 s_24 1
s_10 s_24 1
s_46 w_7 -1
s_33 s_46 1
s_7 s_46 1
s_52 w_6 -1
s_39 s_52 -1
s_22 s_39 -1
s_36 s_39 -1
s_26 s_36 -1
s_14 s_26 1
s_18 s_26 1
s_32 s_36 -1
s_23 s_32 -1
s_27 s_32 -1
s_21 s_52 1
s_43 w_5 -1
s_20 s_43 1
s_41 s_43 -1
s_5 s_41 1
s_38 s_41 -1
s_2 s_38 -1
s_35 s_38 -1
s_1 s_35 -1
s_31 s_35 -1
s_13 s_31 -1
s_26 s_31 -1
s_45 w_4 -1
s_44 s_45 -1
s_19 s_44 -1
s_42 s_44 -1
s_31 s_42 1
s_41 s_42 1
s_4 s_45 1
s_49 w_3 -1
s_34 s_49 -1
s_18 s_34 -1
s_30 s_34 -1
s_16 s_30 1
s_24 s_30 1
s_12 s_49 1
s_48 w_2 -1
s_38 s_48 -1
s_11 s_48 1
s_50 w_8 -1
s_47 s_50 -1
s_28 s_47 -1
s_20 s_28 -1
s_22 s_28 1
s_8 s_47 -1
s_15 s_50 1
s_51 w_1 -1
s_25 s_51 -1
s_1 s_25 -1
s_10 s_25 -1
s_17 s_51 1
"""

var2atoms = {}
for line in scheme.split("\n"):
    if len(line) < 2:
        continue

    var1, var2, sgn = line.split(" ")
    sp_var1 = sp.Symbol(var1)
    sgn = int(sgn)
    if flip and sp_var1 in [sp.Symbol(f"s_{i}") for i in range(1, num_inputs + 1)]:
        idx = int(var1.split("_")[1])
        sp_var1 = sp.Symbol(f"M_{idx}")
        # i = idx % 3
        # j = idx // 3
        # sp_var1 = sp.Symbol(f"A_{i + 1, j + 1}")
        sgn *= -1
    var2atoms[var2] = var2atoms.get(var2, []) + [sgn * sp_var1]

for var in sorted(var2atoms.keys()):
    atoms = var2atoms[var]
    eq = sp.Eq(sp.Symbol(var), sum(atoms))
    print(sp.latex(eq), " \\\\")
