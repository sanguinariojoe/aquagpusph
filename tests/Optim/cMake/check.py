import json

N = 252000

prev_end = 0
with open('ihoc.data', 'r') as f:
    data_txt = f.readlines()
    for i, line in enumerate(data_txt):
        line = line.strip()
        if line == '':
            continue
        start, end = [int(field) for field in line.split(' ')]
        if start == N:
            assert end == N, \
                   f"Invalid ihoc[{i}] = {start}, {end}"
        else:
            assert start == prev_end and end > start, \
                   f"Invalid ihoc[{i}] = {start}, {end}"
            prev_end = end


with open('n_neighs.data', 'r') as f:
    data_txt = f.readlines()
    for i, line in enumerate(data_txt):
        line = line.strip()
        if line == '':
            continue
        legacy, optim = [int(field) for field in line.split(' ')]
        assert legacy == optim, \
               f"Invalid n_neighs[{i}] = {legacy} vs. {optim}"


def get_gputime(fpath, name):
    with open(fpath, 'r') as file:
        data = json.load(file)
    for sample in data['snapshot']['samples']:
        if sample['name'] == name + "::Kernel":
            return sample['end'] - sample['start']
    return float('nan')

legacy = get_gputime('performance.json', "Legacy")
optim = get_gputime('performance.json', "Optim")

assert legacy > optim, \
       f"Optimized version ran {optim - legacy} microsecs slower" + \
       f"({100 * (optim - legacy) / legacy}%)"

print(f"Optimized version ran {legacy - optim} microsecs faster" + \
      f" ({100 * (legacy - optim) / legacy:.1f}%)")
