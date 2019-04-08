FILES = ['set0.00000.dat', 'set1.00000.dat']
N = 5000

count = 0
for fname in FILES:
    with open(fname, 'r') as f:
        data_txt = f.readlines()[18:]
    for line in data_txt:
        if line.strip() == '':
            continue
        assert(int(line.split(',')[1]) == count)
        count += 1

assert(count == 2 * N)
