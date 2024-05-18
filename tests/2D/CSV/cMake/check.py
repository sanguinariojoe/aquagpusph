import filecmp

assert(filecmp.cmp('set0.00000.csv', 'ref.csv', shallow=False))
