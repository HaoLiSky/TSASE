#!/usr/bin/env python

import os

filename = os.path.join(os.path.dirname(__file__), "xyz.glade")
f = open(filename, 'r')
lines = f.readlines()
f.close()
f = open(filename, 'w')
for line in lines:
    if "primary_icon_sensitive" in line:
        continue
    if "secondary_icon_sensitive" in line:
        continue
    if "primary_icon_activatable" in line:
        continue
    if "secondary_icon_activatable" in line:
        continue
    if "use_action_appearance" in line:
        continue
    f.write(line)
f.close()

