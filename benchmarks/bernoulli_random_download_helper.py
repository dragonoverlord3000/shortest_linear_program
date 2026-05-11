import itertools, random, os
random.seed(628)
for n, betax40 in itertools.product(range(12, 29), range(1, 40)):
    M = [[int(40 * random.random() < betax40) for j in range(n)] for i in range(n)]
    os.makedirs("./benchmarks/bernoulli/dataset", exist_ok=True)
    with open(f"./benchmarks/bernoulli/dataset/n_{n}_bx40_{betax40}.txt", "w") as f:
        f.write(f"{n} {n}\n" + "\n".join([" ".join([str(v) for v in row]) for row in M]))
