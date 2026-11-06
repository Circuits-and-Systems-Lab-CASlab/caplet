import sys

def mid(p1, p2):
    return [(p1[0]+p2[0])/2.0, (p1[1]+p2[1])/2.0, (p1[2]+p2[2])/2.0]

def fmt(p):
    return f"{p[0]:g} {p[1]:g} {p[2]:g}"

def refine_panel(line):
    parts = line.strip().split()
    if not parts or parts[0] != 'Q':
        return [line]
    
    header = f"{parts[0]} {parts[1]} {parts[2]}"
    
    # Extract the 4 coordinates (x, y, z)
    p1 = [float(x) for x in parts[3:6]]
    p2 = [float(x) for x in parts[6:9]]
    p3 = [float(x) for x in parts[9:12]]
    p4 = [float(x) for x in parts[12:15]]
    
    # Extract trailing dielectric constants if they exist
    tail = " " + " ".join(parts[15:]) if len(parts) > 15 else ""
    
    # Calculate midpoints and center
    p12 = mid(p1, p2)
    p23 = mid(p2, p3)
    p34 = mid(p3, p4)
    p41 = mid(p4, p1)
    pc = mid(p1, p3)
    
    # Generate the 4 new sub-panels maintaining standard winding order
    sub_panels = [
        f"{header} {fmt(p1)} {fmt(p12)} {fmt(pc)} {fmt(p41)}{tail}",
        f"{header} {fmt(p12)} {fmt(p2)} {fmt(p23)} {fmt(pc)}{tail}",
        f"{header} {fmt(pc)} {fmt(p23)} {fmt(p3)} {fmt(p34)}{tail}",
        f"{header} {fmt(p41)} {fmt(pc)} {fmt(p34)} {fmt(p4)}{tail}"
    ]
    return sub_panels

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python refine_mesh.py <input.qui>")
        sys.exit(1)
        
    with open(sys.argv[1], 'r') as f:
        for line in f:
            for sub in refine_panel(line):
                print(sub)
