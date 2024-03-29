def are_equal(v1, v2, tol=1e-6):
    if len(v1) != len(v2):
        return False
    eps = max([abs(v1[i] - v2[i]) for i in range(len(v1))])
    return eps <= tol


with open("out.00000.dat", "r") as f, \
     open("out.0.00000.dat", "r") as f0, \
     open("out.1.00000.dat", "r") as f1:
    lines = f.readlines()[18:]
    lines0 = f0.readlines()[18:]
    lines1 = f1.readlines()[18:]
    n, n0, n1 = 0, 0, 0
    for i, line in enumerate(lines):
        line = line.replace(',', ' ').strip()
        if line == "":
            continue
        n += 1
        fields = line.split(' ')
        vals = [float(v) for v in fields[:-1]]
        vals.append(int(fields[-1]))
        x = vals[0]
        if x < 0.0:
            line_ref = lines0[n0]
            n0 += 1
        else:
            line_ref = lines1[n1]
            n1 += 1
        line_ref = line_ref.replace(',', ' ').strip()
        fields = line_ref.split(' ')
        vals_ref = [float(v) for v in fields[:-1]]
        vals_ref.append(int(fields[-1]))
        if not are_equal(vals, vals_ref):
            print(f"Particle {i} mismatching:\n\t{vals}\n\t{vals_ref}")
            raise AssertionError
