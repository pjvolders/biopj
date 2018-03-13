class GenomicElement(object):
    def __init__(self, chrom, start, end, strand):
        self.chrom = chrom
        self.start = int(start)  # 1-based start
        self.end = int(end)  # 1-based end
        self.strand = strand

    def __cmp__(self, other):
        if self.end < other.start:
            return -1
        elif self.start > other.end:
            return 1
        # self.end >= other.start && self.start <= other.end
        return 0


class Transcript(GenomicElement):
    def __init__(self, chrom, start, end, strand, name, gene_name = ''):
        """

        :param chrom: str
        :param start: int
        :param end: int
        :param strand: str '+' or '-'
        :param name: str
        :param gene_name: str
        """
        super(Transcript, self).__init__(chrom, start, end, strand)

        self.exons = []

        self.name = name
        self.gene_name = gene_name

        self.cluster = ''

    @property
    def size(self):
        if len(self.exons) > 1:
            my_size = 0
            for e in self.exons:
                my_size += e.size
            return my_size
        return None

    def __str__(self):
        r_val = ''
        for e in self.exons:
            r_val += "\t%s:%d-%d\n"%(e.chrom, e.start, e.end)
        return "%s %s:%d-%d\n"%(self.name, self.chrom, self.start, self.end) + r_val


class Exon(GenomicElement):
    def __init__(self, chrom, start, end, strand, transcript):
        """

        :param chrom: str
        :param start: int
        :param end: int
        :param strand: str
        :param transcript: Transcript
        """
        super(Exon, self).__init__(chrom, start, end, strand)

        self.transcript = transcript

    @property
    def size(self):
        return self.end - self.start + 1

    def __str__(self):
        return self.start + ' - ' + self.end


class GenomicElementList:
    """docstring for GenomicElementList"""

    def __init__(self):
        self.tree = {}  # empty dict

    def add(self, genom_elem):
        if genom_elem.strand not in self.tree:
            self.tree[genom_elem.strand] = {}

        if genom_elem.chrom not in self.tree[genom_elem.strand]:
            self.tree[genom_elem.strand][genom_elem.chrom] = []

        self.tree[genom_elem.strand][genom_elem.chrom].append(genom_elem)

    def sort(self):
        for strand in self.tree:
            for chrom in self.tree[strand]:
                self.tree[strand][chrom].sort(key=lambda ge: ge.start)

    def filter_overlapping(self, genom_elem):
        all_genom_elems = self.tree[genom_elem.strand][genom_elem.chrom]
        genom_elems = []

        for ge in all_genom_elems:
            if ge >= genom_elem:
                if ge == genom_elem:
                    genom_elems.append(ge)
                else:
                    break

        return genom_elems

    def __str__(self):
        r_val = ''
        for strand in self.tree:
            r_val += strand + '\n'
            for chrom in self.tree[strand]:
                r_val += '\t' + chrom + '\n'
                for genom_elem in self.tree[strand][chrom]:
                    r_val += '\t\t' + str(genom_elem) + '\n'

        return r_val
