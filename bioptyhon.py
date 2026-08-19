class Align:
    miss_score = -1
    match_score = +1
    gap_score = -2
    mode = "loacl"

    def __init__(self):
        pass

    def align(self,seq1,seq2):
        if self.mode == "local":
            #smith
        else: