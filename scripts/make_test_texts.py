#!/usr/bin/env python3
"""Generate the three 8 MB test texts used in the paper into data/:
headers.txt (concatenated C headers from /usr/include), dna-rand.txt, dna-rep.txt.
Note: headers.txt depends on the local /usr/include contents."""
import glob, os, random

out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
os.makedirs(out, exist_ok=True)
N = 8_000_000

buf = bytearray()
for f in sorted(glob.glob("/usr/include/**/*.h", recursive=True)):
    try:
        buf += open(f, "rb").read()
    except OSError:
        pass
    if len(buf) >= N:
        break
buf = buf[:N].replace(b"\0", b"\1")
open(os.path.join(out, "headers.txt"), "wb").write(buf)

random.seed(1)
open(os.path.join(out, "dna-rand.txt"), "w").write("".join(random.choice("ACGT") for _ in range(N)))
base = "".join(random.choice("ACGT") for _ in range(100_000))
parts = []
for _ in range(80):
    b = list(base)
    for _ in range(100):
        b[random.randrange(len(b))] = random.choice("ACGT")
    parts.append("".join(b))
open(os.path.join(out, "dna-rep.txt"), "w").write("".join(parts))
print("written to", out)
