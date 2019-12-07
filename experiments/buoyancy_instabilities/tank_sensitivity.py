""" run a series of tank experiments with different viscosities """

import os

exps = [
    ['param.IO["expname"] = "khi_novisc"',
     'param.physics["diff_coef"] = {}'],
    #
    ['param.IO["expname"] = "khi_visc_0.1"',
     'param.physics["diff_coef"] = {"u": 1e-1}']
]

script = "tank.py"

with open(script, "r") as fid:
    lines = fid.readlines()

params = [p.split('=')[0] for p in exps[0]]

for e in exps:
    print('*'*80)
    with open(script, "w") as fid:
        for l in lines:
            param = l.split('=')[0]
            if param in params:
                p = [q for q in e if param in q][0]
                print(p)
                l = p+'\n'

            fid.write(l)

    os.system("python "+script)
