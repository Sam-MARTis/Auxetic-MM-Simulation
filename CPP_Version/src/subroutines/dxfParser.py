import ezdxf
from pathlib import Path
import sys

if __name__ == "__main__":
    dxf_path = Path(sys.argv[1])
    doc = ezdxf.readfile(dxf_path)
    msp = doc.modelspace()

    boundary = []
    cuts = []
    for line in msp.query('LINE[layer=="BOUNDARY"]'):
        boundary.append((line.dxf.start, line.dxf.end))

    for line in msp.query('LINE[layer=="CUTS"]'):
        cuts.append((line.dxf.start, line.dxf.end))

    print(len(boundary))
    for start, end in boundary:
        print(start[0], start[1], end[0], end[1])

    print(len(cuts))
    for start, end in cuts:
        print(start[0], start[1], end[0], end[1])
