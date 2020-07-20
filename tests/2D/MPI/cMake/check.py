PAIRS = [{'in':'particles_0.dat', 'out':'out_1.00000.dat'},
         {'in':'particles_1.dat', 'out':'out_0.00000.dat'},]


def is_equal(a, b, tol=1E-5):
    return abs(a - b) < tol


for pair in PAIRS:
    r = []
    with open(pair['out'], 'r') as f:
        data_txt = f.readlines()[18:]
    for line in data_txt:
        line = line.strip()
        if line == '':
            continue
        line = line.replace(',', '')
        r.append([float(field) for field in line.split(' ')])

    with open(pair['in'], 'r') as f:
        data_txt = f.readlines()
    for i, line in enumerate(data_txt):
        line = line.strip()
        if line == '':
            continue

        for j, field in enumerate(line.split(' ')[:-1]):
            assert(is_equal(r[i][j], float(field)))
