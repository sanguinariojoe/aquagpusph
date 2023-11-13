with open('set0.00000.dat', 'r') as f:
    data_txt = f.readlines()[18:]
    f, f_orig, perms, inv_perms = [], [], [], []
    for i, line in enumerate(data_txt):
        line = line.replace(",", " ").strip()
        if line == '':
            continue
        fields = [int(field) for field in line.split(' ')]
        f.append(fields[0])
        f_orig.append(fields[1])
        perms.append(fields[2])
        inv_perms.append(fields[3])

assert(f == sorted(f))
assert(f == sorted(f_orig))
assert(f_orig != sorted(f_orig))
for i in range(len(f)):
    assert(f[i] == f_orig[perms[i]])
    assert(f[inv_perms[i]] == f_orig[i])
