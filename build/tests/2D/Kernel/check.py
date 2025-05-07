def is_approx(v1, v2, tol=1.001e-5):
    return abs(v1 - v2) < tol

with open('set0.00000.dat', 'r') as f:
    data_txt = f.readlines()[18:]
    for line in data_txt:
        line = line.replace(",", " ").strip()
        if line == '':
            continue
        fields = [float(field) for field in line.split(' ')]
        n = len(fields) // 2
        for i in range(n):
            assert(is_approx(2 * fields[i], fields[i + n]))
