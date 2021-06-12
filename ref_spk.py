import json

class Ref_spk():
    def __init__(self):
        None
        
    def set_path(self, path):
        self.path = path
        
    def read_spk(self):
        with open(self.path, "r") as f:
            data = f.read()
            self.trans = json.loads(data)
        
    def fill_isotope_trans(self, isotope):
        self.transitions = []
        
        try:
            for line in self.trans[isotope].keys():
                en = float(line)
                I = self.trans[isotope][line]
                self.transitions.append([en, I])
        except KeyError:
            print("Error: no entry for {} isotope".format(isotope))
            raise Exception("KeyError")
