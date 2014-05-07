
from multiprocessing import Queue, Process

class xyz_process():
    def __init__(self, title="xyz"):
        self.title = title
        self.qin = Queue()
        self.qout = Queue()
        self.process = Process(target=self.target)
        self.process.daemon = True
        self.process.start()
    def target(self):
        from xyz import xyz, gtk
        _xyz = xyz(self.qin, self.qout, title=self.title)
        gtk.main()
    def put(self, atoms):
        self.qout.put(atoms, False)
    def get(self):
        atoms = None
        while not self.qin.empty():
            atoms = self.qin.get(False)
        return atoms


