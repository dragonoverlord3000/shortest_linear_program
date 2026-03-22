fp = "./schemes-tab/schemes/14x2+55x+4y6+2y4+17y3+2z28+z27/i1w213c48ah-000.tab"
type = "U"
gf2 = True

m = 9 if type == "W" else 23
n = 23 if type == "W" else 9

G = [[0] * n for _ in range(m)]

simplifier = lambda x: abs(x) if gf2 else x

def insert_UV(row, main_idx):
    row_idx = main_idx%3
    for col_idx in range(3):
        G[main_idx//3][3*row_idx + col_idx] = simplifier(row[col_idx])

def insert_W(row, main_idx):
    row_idx = main_idx%3
    for col_idx in range(3):
        G[3*row_idx + col_idx][main_idx//3] = simplifier(row[col_idx])

with open(fp, "r") as f:
    data = f.read().split("\n")

    main_idx = 0
    for line in data:
        if "--" in line or not line:
            continue

        inner_data = line.split("|")
        row = inner_data[{"U": 0, "V": 1, "W": 2}[type]]
        row = row.strip(" ")
        row = list(map(int, row.split()))

        if type in "UV":
            insert_UV(row, main_idx)
        else:
            insert_W(row, main_idx)

        main_idx += 1

f.close()


print("{", end="")
for row_idx in range(m):
    row = G[row_idx]
    print("{", end="")
    for idx in range(len(row)):
        print(row[idx], end="")
        if idx != n - 1:
            print(",", end="")
    print(["},", "}"][row_idx == m-1], end="")
print("}")
