with open("log.proc0.00000.html", 'r') as f:
    txt = f.read()
assert(r'[INFO] (InstallableDemo::_execute): Executing Installed demo tool!' in txt)
assert(r'[INFO] (main): Simulation finished OK (t = 0 s)' in txt)
