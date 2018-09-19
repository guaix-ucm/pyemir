
import glob
import itertools
import astropy.io.fits as fits
import yaml


counter_k1 = ['OBSBLOCK', 'IMGOBBL', 'EXP', 'FRSEC']
counter_k2 = ['NOBSBLCK', 'NIMGOBBL', 'NEXP', 'NFRSEC']
target_k = ['DATE-OBS', 'OBSMODE', 'OBJECT', 'OBSTYPE']
counters_v1 = ['EXP', 'NEXP']


class ObsBlockNode(object):
    def __init__(self, maxlen=1):
        self.id = 0
        self.nodes = []
        self.maxlen = maxlen
        self.parent = None
        self.mode = 'PARENT'
        self.index = 0
        self.frames = []
        self.value = 0
        self.label = ""

    def __len__(self):
        return len([x for x in self.nodes if x is not None])

    def capacity(self):
        return sum(node.capacity() for node in self.nodes)

    def occupancy(self):
        return sum(node.occupancy() for node in self.nodes)


class ImageNode(ObsBlockNode):
    def __init__(self, maxlen=1):
        super(ImageNode, self).__init__(maxlen)
        self.frames = None
        self.mode = "ImageNode"

    def capacity(self):
        return 1

    def occupancy(self):
        return int(self.frames is not None)


def create_branch(a, b, c, d=0):
    obid = 1
    nodea = ObsBlockNode(maxlen=a)
    nodea.label = "OO"
    nodea.index = 1
    nodea.id = obid
    obid += 1
    for idx1 in range(a):
        nodeb = ObsBlockNode(maxlen=b)
        nodeb.parent = nodea
        nodeb.index = idx1 + 1
        nodeb.label = "L2"
        nodeb.id = obid
        obid +=1
        nodea.nodes.append(nodeb)
        for idx2 in range(b):
            nodec = ObsBlockNode(maxlen=c)
            nodec.label = "L3"
            nodec.index = idx2 + 1
            nodec.parent = nodeb
            nodec.id = obid
            obid += 1
            nodeb.nodes.append(nodec)
            for idx3 in range(c):
                noded = ImageNode(maxlen=d)
                noded.label = "Image"
                noded.index = idx3 + 1
                noded.parent = nodec
                noded.id = obid
                obid += 1
                nodec.nodes.append(noded)
                for idx4 in range(d):
                    nodee = ImageNode(maxlen=1)
                    nodee.index = idx4 + 1
                    nodee.label = "ImageRaw"
                    nodee.parent = noded
                    nodee.id = obid
                    obid += 1
                    noded.nodes.append(nodee)
    return nodea


def find_node(tree, image, counters):
    if counters:
        counter, rest = counters[0], counters[1:]
        # counter, *rest = counters
        # select node
        idx = counter - 1
        subtree = tree.nodes[idx]
        return find_node(subtree, image, rest)
    return tree


def visit_tree_p(tree, h='-'):

    if(isinstance(tree, ImageNode)):
        return
    for node in tree.nodes:
        visit_tree_p(node, h='-'+h)
    print(h+'>', 'id=', tree.id, 'index=', tree.index, 'mode', tree.mode, 'capacity=',tree.maxlen)


def visit_tree_r(tree, prev):

    res = {'id': tree.id, 'children': [], 'instrument': 'EMIR', 'mode': tree.mode,
           'frames': tree.frames}
    for node in tree.nodes:
        if isinstance(node, ImageNode):
            continue
        subres = visit_tree_r(node, prev)
        res['children'].append(subres['id'])
    prev.append(res)
    return res


if __name__ == '__main__':
    import pickle

    oric = False

    if oric:
        values = {}
        for fname in glob.glob('*EMIR*.fits'):
            hdr = fits.getheader(fname)
            ivalues = {}
            for k in itertools.chain(counter_k1, counter_k2, target_k):
                ivalues[k] = hdr[k]
            values[fname] = ivalues
        with open('inter.pkl', 'wb') as fd:
            pickle.dump(values, fd)
    else:
        with open('inter.pkl', 'rb') as fd:
            values = pickle.load(fd)

    # print(values)
    trees = []
    tree = create_branch(0, 0, 0)
    images_s = sorted(values.keys())
    image_si = iter(images_s)
    image = next(image_si)
    counter = 1
    errorc = 0
    create_new = False

    while True:
        v = values[image]
        this = [v[m] for m in counter_k2]
        that = [v[m] for m in counter_k1]
        other = [v[m] for m in target_k]
        #
        # print(image, this[:3], that[:3], other)
        # print(image, other)

        if create_new or (tree.occupancy() >= tree.capacity()):
            create_new = False
            print('block with counters', image, this[:3])
            ntotal = this[0] * this[1] * this[2]
            print('create new tree', 'expecting', ntotal, 'image(s)')
            tree = create_branch(*this[:3])
            tree.label = counter
            counter += 1
            trees.append(tree)
        # print('capacity', tree.capacity())
        # print('occ', tree.occupancy())
        # room = tree.capacity() - tree.occupancy()
        # print('room', room)
        print('add image', image, that[:3])
        try:
            m = find_node(tree, image, that[:3])
            if m.frames is None:
                m.frames = [image]
                m.parent.mode = values[image]['OBSMODE']
                m.parent.frames.append(image)
            else:
                create_new = True
                raise IndexError('node already used')
            # go to next
            try:
                image = next(image_si)
            except StopIteration:
                break
        except IndexError as error:
            print(error)
            print('cannot match image, try next loop')
            errorc += 1
            if errorc > 10:
                raise

    prev = []
    for tree in trees:
        visit_tree_r(tree, prev)

    with open('obs.yaml', 'w') as fd:
        yaml.dump_all(prev, fd)
    #visit_tree_p(tree)