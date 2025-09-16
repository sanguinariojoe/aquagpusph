with open('log.proc0.00000.html', 'r') as f:
    assert("invalid document structure" in f.read())
with open('log.proc0.00001.html', 'r') as f:
    assert("expected end of tag 'UnclosedDOM'" in f.read())
with open('log.proc0.00002.html', 'r') as f:
    assert("'<' character cannot be used in attribute" in f.read())
