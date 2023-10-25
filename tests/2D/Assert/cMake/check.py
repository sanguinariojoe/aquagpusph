with open("log.proc0.00000.html", 'r') as f:
    txt = f.read()
assert(r'[ERROR] (CalcServer::exec_status_check): Tool "Bad assertion" exited with the error code -1.' in txt)
assert(r'[WARNING] (CalcServer::solver): Skipping "skipped tool" due to dependency errors.' in txt)
assert(r'[INFO] (main): Simulation finished abnormally (t = 0 s)' in txt)
