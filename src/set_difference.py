import sys


class segment:
    def __init__(self, input_str):
        self.input_str = input_str
        lst = input_str.strip().split()
        self.start = int(float(lst[5]))
        self.end = int(float(lst[6]))
        self.haplotype1 = lst[1]
        self.haplotype2 = lst[3]
        self.tsk_id1 = lst[0]
        self.tsk_id2 = lst[2]
    
    def __eq__(self, other):
        if (self.tsk_id1 == other.tsk_id1) and (self.tsk_id2 == other.tsk_id2):
            return True
        else:
            return False
        

def main():
    hi_input = sys.argv[1]
    hi_ps_input = sys.argv[2]
    f = open(hi_input)
    hi_segments = [segment(i) for i in f.readlines()]
    f.close()
    f = open(hi_ps_input)
    hi_ps_segments = [segment(i) for i in f.readlines()]
    f.close()
    f = open("set_difference.txt", 'w')
    for seg1 in hi_ps_segments:
        found = False
        for seg2 in hi_segments:
            if seg1 == seg2:
                found = True
        if not found:
            f.write(seg1.input_str)
    f.close()



if __name__ == "__main__":
    main()






