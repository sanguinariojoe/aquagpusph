v_org = []
with open('particles.dat', 'r') as f:
    data_txt = f.readlines()
    for line in data_txt:
        line = line.replace(",", " ").strip()
        if line == '':
            continue
        fields = [float(field) for field in line.split(' ')]
        v_org.append(int(fields[2]))

def accumu(v):
    total = 0
    yield total
    for x in v[:-1]:
        total += x
        yield total
v_org = list(accumu(v_org))

v_dst = []
with open('set0.00000.dat', 'r') as f:
    data_txt = f.readlines()[18:]
    for line in data_txt:
        line = line.replace(",", " ").strip()
        if line == '':
            continue
        fields = [float(field) for field in line.split(' ')]
        v_dst.append(int(fields[2]))

for i in range(len(v_org)):
    assert(v_org[i] == v_dst[i])
