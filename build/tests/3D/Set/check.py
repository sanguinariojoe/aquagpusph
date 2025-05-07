values = (0, 0, 0, 0, 0.2)

with open('set0.00000.dat', 'r') as f:
    data_txt = f.readlines()[18:]
    for line in data_txt:
        line = line.replace(",", " ").strip()
        if line == '':
            continue
        fields = [float(field) for field in line.split(' ')]
        for i in range(len(fields)):
            assert(fields[i] == values[i])
