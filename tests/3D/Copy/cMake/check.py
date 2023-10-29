with open('set0.00000.dat', 'r') as f:
    data_txt = f.readlines()[18:]
    for line in data_txt:
        line = line.replace(",", " ").strip()
        if line == '':
            continue
        fields = [float(field) for field in line.split(' ')]
        n = len(fields) // 2
        for i in range(n):
            assert(fields[i] == fields[i + n])
