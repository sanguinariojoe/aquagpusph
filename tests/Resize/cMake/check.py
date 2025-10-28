with open('pre.out', 'r') as fpre:
    assert(len(fpre.readlines()) == 128)

with open('post.out', 'r') as fpost:
    assert(len(fpost.readlines()) == 500)

with open('shrink.out', 'r') as fshrink:
    assert(len(fshrink.readlines()) == 500)
