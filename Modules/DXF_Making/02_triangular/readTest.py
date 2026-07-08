import ezdxf
from pathlib import Path
import sys
if __name__ == "__main__":
    dxf_path = Path(sys.argv[1])
    doc = ezdxf.readfile(dxf_path)
    msp = doc.modelspace()

    boundary = []
    numBoundaryPoints = 0
    cuts = []
    numCutPairs = 0
    for line in msp.query('LINE[layer=="BOUNDARY"]'):
        boundary.append((line.dxf.start, line.dxf.end))
        numBoundaryPoints += 1

    for line in msp.query('LINE[layer=="CUTS"]'):
        cuts.append((line.dxf.start, line.dxf.end))
        numCutPairs += 1

    # print("Boundary lines:")
    print(numBoundaryPoints)
    for start, end in boundary:
        print(start[0], start[1])
    
    
    # print("\nCut lines:")
    print(numCutPairs)
    for start, end in cuts:
        print(start[0], start[1],end[0], end[1])
    
