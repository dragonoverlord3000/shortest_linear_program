import itertools, random, os
random.seed(628)
for n, betax10 in itertools.product(range(15, 21), range(1, 10)):
    M = [[0] * n for _ in range(n)]
    for i, j in itertools.product(range(n), range(n)):
        M[i][j] = int(10 * random.random() < betax10)
    os.makedirs("./benchmarks/bernoulli/dataset", exist_ok=True)
    with open(f"./benchmarks/bernoulli/dataset/n_{n}_bx10_{betax10}.txt", "w") as f:
        f.write(f"{n} {n}\n" + "\n".join([" ".join([str(v) for v in row]) for row in M]))
