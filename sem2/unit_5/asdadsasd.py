



unc = {
    "Earth" :  {"period" :[], "perihelion":[], "aphelion":[]},
    "Moon"  :  {"period" :[], "perihelion":[], "aphelion":[]},
    "Mercury": {"period" :[], "perihelion":[], "aphelion":[]}
}


unc["Earth"]['period'].append(1)
print(unc)

unc["Earth"]['period'].append(2)
print(f"\n\n{unc}")