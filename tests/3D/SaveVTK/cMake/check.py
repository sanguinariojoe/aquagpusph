import filecmp

assert(filecmp.cmp('set0.00000.vtu', 'ref.vtu', shallow=False))
