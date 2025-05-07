def read(fpath):
    data = []
    with open(fpath, "r") as f:
        lines = f.readlines()
        for l in lines:
            l = l.strip()
            comment = l.find('#')
            if comment != -1:
                l = l[:comment]
            if l == "":
                continue
            l = l.replace(",", " ")
            while l.find("  ") != -1:
                l = l.replace("  ", " ")
            fields = l.split(" ")
            data.append([float(field) for field in fields])
    return data


def are_equal(x, y, tol=1e-5):
    return abs(x - y) < tol

org = read("Fluid.dat")
predump = read("pre_predictor_dump.00000.txt")
postdump = read("post_predictor_dump.00000.txt")
output = read("output.00000.dat")

for i in range(len(predump)):
    for j in range(len(predump[i])):
        assert are_equal(org[i][j], predump[i][j])
        assert are_equal(org[i][j], postdump[i][j])
        assert are_equal(org[i][j], output[i][j])
