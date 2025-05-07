with open('set0.00000.vtu', mode='r') as fa, open('ref.vtu', mode='r') as fb:
    assert(fa.read() == fb.read())
