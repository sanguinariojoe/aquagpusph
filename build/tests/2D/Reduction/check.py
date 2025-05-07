def is_approx(v1, v2, tol=1e-5):
    return abs(v1 - v2) < tol

r_sum = [0.0, 0.0]
r_max = [-float("inf"), -float("inf")]
with open('set0.00000.dat', 'r') as f:
    data_txt = f.readlines()[18:]
    for line in data_txt:
        line = line.replace(",", " ").strip()
        if line == '':
            continue
        fields = [float(field) for field in line.split(' ')]
        for i, field in enumerate(fields):
            r_sum[i] += field
            r_max[i] = max(r_max[i], field)

with open('vars.out', 'r') as f:
    line = f.readlines()[1]
    line = line.replace(",", " ").replace("(", "").replace(")", "").strip()
    fields = [float(field) for field in line.split(' ')]
    n_fields = len(fields) // 2
    for i in range(n_fields):
        assert(is_approx(r_sum[i], fields[i], tol=1e-2))
        assert(is_approx(r_max[i], fields[i + n_fields]))
