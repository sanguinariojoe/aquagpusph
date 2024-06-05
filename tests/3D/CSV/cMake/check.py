def is_equal(a, b, tol=1E-7):
    return abs(a - b) < tol


with open('set0.00000.csv', 'r') as f1, open('ref.csv', 'r') as f2:
    for l1, l2 in zip(f1.readlines(), f2.readlines()):
        l1 = l1.strip()
        l2 = l2.strip()
        if l1 == l2:
            # We cannot dream with a better match
            continue
        fields1 = [float(field) for field in l1.split(",")]
        fields2 = [float(field) for field in l2.split(",")]
        for field1, field2 in zip(fields1, fields2):
            assert(is_equal(field1, field2))
