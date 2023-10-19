with open("vars.out", 'r') as f:
    data_txt = f.readlines()[1]
assert(data_txt.strip() == '-1 -4 (-16,-0.75,5,6)')
