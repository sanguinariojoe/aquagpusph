with open("log.proc0.00000.html", 'r') as f:
    txt = f.read()
assert(r'[ERROR] (Assert::_solve): Assertion error. The expression "i < 100" on tool "Stop" is false' in txt)
assert(r'[INFO] (main): Simulation finished abnormally (t = 0 s)' in txt)
